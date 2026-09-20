"""Study-family evidence classification and synthesis."""

from __future__ import annotations

import numpy as np

import pandas as pd

EXPECTED_CONTEXTS = {
    "FAMILY_1": ["HIEC6", "HaCaT"],
    "FAMILY_2": ["iPSC-Hep", "iPSC-CM"],
    "FAMILY_3": ["liver_spheroid_24h", "liver_spheroid_240h"],
}

HETEROGENEITY_LABELS = {
    "FAMILY_1": "STUDY_CONTEXT_CELL_LINE_BATCH",
    "FAMILY_2": "IPSC_CONTEXT_PROTOCOL_DURATION",
    "FAMILY_3": "DURATION_BATCH_CONTEXT",
}

def to_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.upper().isin(["TRUE", "T", "1"])

def classify_evidence(frame: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    t = cfg["thresholds"]
    out = frame.copy()
    for col in [
        "camera_fdr",
        "score_fdr",
        "score_standardized_effect",
        "camera_p_value",
        "score_p_value",
    ]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out["eligible"] = to_bool(out["eligible"])
    out["camera_direction"] = out["camera_direction"].astype(str).str.upper()
    out["score_direction"] = out["score_direction"].astype(str).str.upper()
    direction_agreement = out["camera_direction"].eq(out["score_direction"])
    finite_score = np.isfinite(out["score_standardized_effect"])
    out["strong"] = (
        out["eligible"]
        & (out["camera_fdr"] <= t["camera_strong_fdr"])
        & finite_score
        & (out["score_standardized_effect"].abs() >= t["score_abs_standardized_effect"])
        & (out["score_fdr"] <= t["score_fdr"])
        & direction_agreement
    )
    camera_weak = out["camera_fdr"] <= t["pathway_weak_fdr"]
    score_weak = out["score_fdr"] <= t["score_fdr"]
    method_conflict = (
        out["eligible"]
        & camera_weak
        & score_weak
        & ~direction_agreement
    )
    out["weak"] = out["eligible"] & ~out["strong"] & ~method_conflict & (camera_weak | score_weak)
    out["evidence_direction"] = ""
    out.loc[out["strong"], "evidence_direction"] = out.loc[out["strong"], "camera_direction"]
    use_camera = out["weak"] & camera_weak
    out.loc[use_camera, "evidence_direction"] = out.loc[use_camera, "camera_direction"]
    use_score = out["weak"] & ~camera_weak & score_weak
    out.loc[use_score, "evidence_direction"] = out.loc[use_score, "score_direction"]
    out["evidence_class"] = "NO_EVIDENCE"
    out.loc[~out["eligible"], "evidence_class"] = "NOT_ASSESSED"
    out.loc[method_conflict, "evidence_class"] = "METHOD_CONFLICT"
    out.loc[out["weak"], "evidence_class"] = "WEAK"
    out.loc[out["strong"], "evidence_class"] = "STRONG"
    return out

def build_heterogeneity(all_results: pd.DataFrame, context: pd.DataFrame) -> pd.DataFrame:
    interactions = all_results[all_results["contrast_type"] == "HETEROGENEITY"].copy()
    rows: list[dict] = []
    pathways = sorted(context["pathway"].unique())
    for family, expected_contexts in EXPECTED_CONTEXTS.items():
        for pathway in pathways:
            c = context[(context["family_id"] == family) & (context["pathway"] == pathway)]
            i = interactions[(interactions["family_id"] == family) & (interactions["pathway"] == pathway)]
            formal = bool(i["strong"].any())
            primary_strong = c[c["strong"]]
            directions = set(primary_strong["evidence_direction"])
            opposing = "UP" in directions and "DOWN" in directions
            supported = (formal and not primary_strong.empty) or opposing
            basis = ""
            if formal and not primary_strong.empty:
                basis = "FORMAL_INTERACTION"
            if opposing:
                basis = "OPPOSING_STRONG_CONTEXTS" if not basis else basis + "+OPPOSING_STRONG_CONTEXTS"
            rows.append(
                {
                    "family_id": family,
                    "pathway": pathway,
                    "expected_context_count": len(expected_contexts),
                    "formal_interaction_strong": formal,
                    "opposing_strong_contexts": opposing,
                    "supported_context_heterogeneity": supported,
                    "basis": basis,
                    "interpretation_label": HETEROGENEITY_LABELS[family],
                    "interaction_context": "|".join(i["context"].astype(str).unique()),
                    "interaction_camera_fdr_min": pd.to_numeric(i["camera_fdr"], errors="coerce").min(),
                    "interaction_score_fdr_min": pd.to_numeric(i["score_fdr"], errors="coerce").min(),
                    "interaction_score_abs_std_max": pd.to_numeric(
                        i["score_standardized_effect"], errors="coerce"
                    ).abs().max(),
                }
            )
    return pd.DataFrame(rows)

def build_family_evidence(context: pd.DataFrame, heterogeneity: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    pathways = sorted(context["pathway"].unique())
    for family, expected_contexts in EXPECTED_CONTEXTS.items():
        for pathway in pathways:
            sub = context[(context["family_id"] == family) & (context["pathway"] == pathway)].copy()
            sub = sub.set_index("context").reindex(expected_contexts).reset_index()
            strong_dirs = set(sub.loc[sub["strong"].fillna(False), "evidence_direction"])
            if "UP" in strong_dirs and "DOWN" in strong_dirs:
                vote = "CONFLICT"
            elif "UP" in strong_dirs:
                vote = "UP"
            elif "DOWN" in strong_dirs:
                vote = "DOWN"
            else:
                vote = ""
            weak_dirs = set(sub.loc[sub["weak"].fillna(False), "evidence_direction"])
            weak_vote = ""
            if not vote and len(weak_dirs) == 1:
                weak_vote = next(iter(weak_dirs))
            same_support = 0
            r1 = False
            if vote in {"UP", "DOWN"}:
                same_support = int(
                    (
                        sub["evidence_class"].isin(["STRONG", "WEAK"])
                        & sub["evidence_direction"].eq(vote)
                    ).sum()
                )
                opposite_strong = bool(
                    (
                        sub["strong"].fillna(False)
                        & sub["evidence_direction"].ne(vote)
                    ).any()
                )
                r1 = same_support >= 2 and not opposite_strong
            hrow = heterogeneity[
                (heterogeneity["family_id"] == family)
                & (heterogeneity["pathway"] == pathway)
            ].iloc[0]
            model_restricted = (
                vote in {"UP", "DOWN"}
                and same_support == 1
                and not bool(hrow["supported_context_heterogeneity"])
            )
            rows.append(
                {
                    "family_id": family,
                    "pathway": pathway,
                    "eligible_contexts": int(sub["eligible"].fillna(False).sum()),
                    "strong_contexts": int(sub["strong"].fillna(False).sum()),
                    "weak_contexts": int(sub["weak"].fillna(False).sum()),
                    "family_vote": vote,
                    "family_weak_support": weak_vote,
                    "R1_within_family_consistency": r1,
                    "model_restricted": model_restricted,
                    "supported_context_heterogeneity": bool(
                        hrow["supported_context_heterogeneity"]
                    ),
                    "context_states": ";".join(
                        f"{row.context}:{row.evidence_class if pd.notna(row.evidence_class) else 'NOT_ASSESSED'}"
                        + (
                            f"_{row.evidence_direction}"
                            if pd.notna(row.evidence_direction) and row.evidence_direction
                            else ""
                        )
                        for row in sub.itertuples()
                    ),
                }
            )
    return pd.DataFrame(rows)

def build_cross_family(family: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for pathway, sub in family.groupby("pathway", sort=True):
        sub = sub.set_index("family_id").reindex(EXPECTED_CONTEXTS).reset_index()
        votes = dict(zip(sub["family_id"], sub["family_vote"].fillna("")))
        weak = dict(zip(sub["family_id"], sub["family_weak_support"].fillna("")))
        up = [f for f, v in votes.items() if v == "UP"]
        down = [f for f, v in votes.items() if v == "DOWN"]
        if len(up) >= 2 and not down:
            direction, contributing = "UP", up
        elif len(down) >= 2 and not up:
            direction, contributing = "DOWN", down
        else:
            direction, contributing = "", []
        r2 = bool(direction)
        third_support = False
        r1_pair = False
        if r2:
            noncontributing = [f for f in EXPECTED_CONTEXTS if f not in contributing]
            third_support = any(
                votes.get(f) == direction or weak.get(f) == direction
                for f in noncontributing
            )
            r1_pair = sum(
                bool(
                    sub.loc[
                        sub["family_id"].eq(f), "R1_within_family_consistency"
                    ].iloc[0]
                )
                for f in contributing
            ) >= 2
        r3 = r2 and (third_support or r1_pair)
        grade = "R3" if r3 else ("R2" if r2 else "NOT_REPRODUCIBLE")
        if up and down:
            grade = "FAMILY_DIRECTION_CONFLICT"
        rows.append(
            {
                "pathway": pathway,
                "reproducibility_grade": grade,
                "reproducible_direction": direction,
                "R2": r2,
                "R3": r3,
                "strong_family_count": len(contributing),
                "contributing_families": "|".join(contributing),
                "third_family_same_direction_support": third_support,
                "two_contributing_families_R1": r1_pair,
                "FAMILY_1_vote": votes.get("FAMILY_1", ""),
                "FAMILY_2_vote": votes.get("FAMILY_2", ""),
                "FAMILY_3_vote": votes.get("FAMILY_3", ""),
                "FAMILY_1_R1": bool(
                    sub.loc[
                        sub["family_id"].eq("FAMILY_1"), "R1_within_family_consistency"
                    ].iloc[0]
                ),
                "FAMILY_2_R1": bool(
                    sub.loc[
                        sub["family_id"].eq("FAMILY_2"), "R1_within_family_consistency"
                    ].iloc[0]
                ),
                "FAMILY_3_R1": bool(
                    sub.loc[
                        sub["family_id"].eq("FAMILY_3"), "R1_within_family_consistency"
                    ].iloc[0]
                ),
            }
        )
    return pd.DataFrame(rows)

def build_direct(all_results: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    direct_types = {
        "DIRECT_COMPARATOR",
        "DIRECT_COMPARATOR_EFFECT_PATTERN_ONLY",
        "CHEMICAL_DOSE_INTERACTION_BASIS",
    }
    direct = all_results[
        all_results["contrast_type"].isin(direct_types)
        & all_results["role"].isin(["PRIMARY"])
        & all_results["family_id"].isin(["FAMILY_1", "FAMILY_2"])
    ].copy()
    direct["direct_support"] = ""
    direct.loc[
        direct["strong"]
        & direct["contrast_type"].isin(
            ["DIRECT_COMPARATOR", "DIRECT_COMPARATOR_EFFECT_PATTERN_ONLY"]
        ),
        "direct_support",
    ] = "SUPPORTED_DIRECT"
    interaction_rows = direct[
        direct["contrast_type"] == "CHEMICAL_DOSE_INTERACTION_BASIS"
    ]
    for idx, row in interaction_rows.iterrows():
        anchored = direct[
            (direct["family_id"] == row["family_id"])
            & (direct["context"] == row["context"])
            & (direct["pathway"] == row["pathway"])
            & (direct["contrast_type"] == "DIRECT_COMPARATOR")
            & direct["strong"]
        ]
        if bool(row["strong"]) and not anchored.empty:
            direct.loc[idx, "direct_support"] = "SUPPORTED_INTERACTION_WITH_DIRECT_DOSE_ANCHOR"
    supported = direct[direct["direct_support"].ne("")].copy()
    return direct, supported
