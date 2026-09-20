"""Synthesize primary common-dose and dose-sensitivity evidence."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import pandas as pd


ARMS = {
    "primary_common_five": {
        "context_type": "PFTeDA_CONTEXT_SUMMARY",
        "interaction_type": "HETEROGENEITY",
        "interpretation": "Primary prespecified equal-weight summary across the five doses shared by 24 h and 240 h.",
    },
    "sensitivity_all_available": {
        "context_type": "PFTeDA_F3_ALL_AVAILABLE_DOSE_SUMMARY",
        "interaction_type": "HETEROGENEITY_SENSITIVITY_FOUR_DOSES",
        "interpretation": "All-available-dose context summaries: six doses at 24 h and five doses at 240 h; duration/batch interaction evaluated across the lowest four shared doses.",
    },
    "sensitivity_lowest_four": {
        "context_type": "PFTeDA_F3_LOWEST_FOUR_SHARED_DOSE_SUMMARY",
        "interaction_type": "HETEROGENEITY_SENSITIVITY_FOUR_DOSES",
        "interpretation": "Equal-weight summary across the lowest four shared doses that use the ctrl01 solvent-control group.",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--evidence-module", required=True)
    return parser.parse_args()


def load_evidence(path: Path):
    spec = importlib.util.spec_from_file_location("evidence_module", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load evidence module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_context(all_results: pd.DataFrame, f3_type: str) -> pd.DataFrame:
    non_f3 = all_results[
        (all_results["role"] == "PRIMARY")
        & (all_results["family_id"] != "FAMILY_3")
        & (all_results["contrast_type"] == "PFTeDA_CONTEXT_SUMMARY")
    ]
    f3 = all_results[
        (all_results["role"] == "PRIMARY")
        & (all_results["family_id"] == "FAMILY_3")
        & (all_results["contrast_type"] == f3_type)
    ]
    context = pd.concat([non_f3, f3], ignore_index=True)
    expected = {
        ("FAMILY_1", "HIEC6"),
        ("FAMILY_1", "HaCaT"),
        ("FAMILY_2", "iPSC-Hep"),
        ("FAMILY_2", "iPSC-CM"),
        ("FAMILY_3", "liver_spheroid_24h"),
        ("FAMILY_3", "liver_spheroid_240h"),
    }
    observed = set(zip(context["family_id"], context["context"]))
    if observed != expected:
        raise RuntimeError(
            f"Context mismatch for {f3_type}: missing={expected-observed}; extra={observed-expected}"
        )
    return context


def build_interaction_source(all_results: pd.DataFrame, f3_interaction_type: str) -> pd.DataFrame:
    base = all_results[
        ~(
            (all_results["family_id"] == "FAMILY_3")
            & (all_results["role"] == "FORMAL_HETEROGENEITY")
        )
    ].copy()
    f3 = all_results[
        (all_results["family_id"] == "FAMILY_3")
        & (all_results["role"] == "FORMAL_HETEROGENEITY")
        & (all_results["contrast_type"] == f3_interaction_type)
    ].copy()
    if f3.empty:
        raise RuntimeError(f"No Family 3 interaction rows for {f3_interaction_type}")
    f3["contrast_type"] = "HETEROGENEITY"
    return pd.concat([base, f3], ignore_index=True)


def main() -> None:
    args = parse_args()
    run_dir = Path(args.run_dir)
    cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    evidence = load_evidence(Path(args.evidence_module))

    raw = pd.read_csv(run_dir / "pathway_results" / "ALL_PATHWAY_RESULTS.csv")
    classified = evidence.classify_evidence(raw, cfg)
    classified.to_csv(
        run_dir / "pathway_results" / "ALL_PATHWAY_RESULTS_WITH_EVIDENCE_CLASS.csv",
        index=False,
    )

    direct, supported_direct = evidence.build_direct(classified)
    context_rows: list[pd.DataFrame] = []
    family_rows: list[pd.DataFrame] = []
    cross_rows: list[pd.DataFrame] = []
    heterogeneity_rows: list[pd.DataFrame] = []
    summaries: dict[str, dict] = {}

    for arm, definition in ARMS.items():
        out_dir = run_dir / arm
        out_dir.mkdir(parents=True, exist_ok=True)
        context = build_context(classified, definition["context_type"])
        interaction_source = build_interaction_source(
            classified, definition["interaction_type"]
        )
        heterogeneity = evidence.build_heterogeneity(interaction_source, context)
        family = evidence.build_family_evidence(context, heterogeneity)
        cross = evidence.build_cross_family(family)

        context.to_csv(out_dir / "CONTEXT_PATHWAY_EVIDENCE.csv", index=False)
        family.to_csv(out_dir / "FAMILY_PATHWAY_EVIDENCE.csv", index=False)
        cross.to_csv(out_dir / "CROSS_FAMILY_REPRODUCIBILITY.csv", index=False)
        cross[cross["R2"]].to_csv(out_dir / "R2_R3_REPRODUCIBLE_MODULES.csv", index=False)
        heterogeneity.to_csv(out_dir / "CONTEXT_HETEROGENEITY_ASSESSMENT.csv", index=False)
        heterogeneity[heterogeneity["supported_context_heterogeneity"]].to_csv(
            out_dir / "SUPPORTED_CONTEXT_HETEROGENEITY.csv", index=False
        )
        direct.to_csv(out_dir / "PFTeDA_PFOA_DIRECT_PROGRAMS_ALL.csv", index=False)
        supported_direct.to_csv(
            out_dir / "SUPPORTED_PFTeDA_PFOA_DIFFERENTIAL_PROGRAMS.csv", index=False
        )

        for frame, collector in [
            (context, context_rows),
            (family, family_rows),
            (cross, cross_rows),
            (heterogeneity, heterogeneity_rows),
        ]:
            tagged = frame.copy()
            tagged.insert(0, "analysis_arm", arm)
            collector.append(tagged)

        reproducible = cross[cross["R2"]][
            ["pathway", "reproducibility_grade", "reproducible_direction"]
        ].to_dict("records")
        emt = cross[cross["pathway"] == "Epithelial Mesenchymal Transition"].iloc[0]
        summaries[arm] = {
            "definition": definition,
            "reproducible_modules": reproducible,
            "emt_grade": emt["reproducibility_grade"],
            "emt_family3_vote": emt["FAMILY_3_vote"],
            "emt_family3_r1": bool(emt["FAMILY_3_R1"]),
            "r2_or_higher_count": int(cross["R2"].sum()),
            "r3_count": int(cross["R3"].sum()),
            "supported_context_heterogeneity_count": int(
                heterogeneity["supported_context_heterogeneity"].sum()
            ),
            "supported_direct_unique_program_count": int(
                supported_direct["pathway"].nunique()
            ),
        }
        (out_dir / "ARM_SUMMARY.json").write_text(
            json.dumps(summaries[arm], indent=2), encoding="utf-8"
        )

    all_context = pd.concat(context_rows, ignore_index=True)
    all_family = pd.concat(family_rows, ignore_index=True)
    all_cross = pd.concat(cross_rows, ignore_index=True)
    all_heterogeneity = pd.concat(heterogeneity_rows, ignore_index=True)
    all_context.to_csv(run_dir / "FAMILY3_CONTEXT_EVIDENCE_SENSITIVITY.csv", index=False)
    all_family.to_csv(run_dir / "FAMILY3_FAMILY_EVIDENCE_SENSITIVITY.csv", index=False)
    all_cross.to_csv(run_dir / "FAMILY3_REPRODUCIBILITY_SENSITIVITY.csv", index=False)
    all_heterogeneity.to_csv(run_dir / "FAMILY3_HETEROGENEITY_SENSITIVITY.csv", index=False)

    comparison = all_cross.pivot_table(
        index="pathway",
        columns="analysis_arm",
        values="reproducibility_grade",
        aggfunc="first",
    ).reset_index()
    comparison.to_csv(run_dir / "FAMILY3_GRADE_COMPARISON.csv", index=False)
    (run_dir / "FAMILY3_SENSITIVITY_SUMMARY.json").write_text(
        json.dumps(summaries, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
