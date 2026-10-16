# Snakemake pipeline. Ingest (real ecosystem rows) -> validate -> modules.
# Run: snakemake --cores 1   (R rules need the R/Bioconductor stack)

rule all:
    input:
        "data/validated/.stamp",
        "results/m1/corrected_or.csv",
        "results/m1/mecor_report.txt",
        "results/m2/dag_adjustment.txt",
        "results/m2/meta_regression.txt",
        "results/m2/baggr_summary.txt",
        "results/m3/occupancy_summary.csv",
        "results/m3/elicited_beta_priors.csv",
        "results/m4/dormancy_grid.csv",
        "results/m4/oncosimulr_diversity.csv",
        "results/m5/cif_brain_first_mut.csv",
        "results/m5/cif_brain_first_wt.csv",
        "results/m6/er_activity.csv"

rule ingest_ecosystem:
    output:
        "data/raw/study_table.parquet",
        "data/raw/meta_anchor.json"
    shell:
        "python data/ingest_ecosystem.py"

rule validate_schema:
    input:
        "data/raw/detection_history.parquet",
        "data/raw/study_table.parquet",
        "data/raw/first_progression.parquet",
        "data/raw/expression_samples.parquet",
        "data/raw/expression.parquet",
        "data/raw/meta_anchor.json"
    output:
        "data/validated/.stamp"
    shell:
        "python data/validate.py"

rule m1_measurement_error:
    input:
        "data/raw/study_table.parquet",
        "data/raw/meta_anchor.json"
    output:
        "results/m1/corrected_or.csv"
    shell:
        "Rscript modules/m1_measurement_error/simex_correction.R"

rule m2_dag:
    output:
        "results/m2/dag_adjustment.txt"
    shell:
        "Rscript modules/m2_dag_meta/dag.R"

rule m2_meta_regression:
    input:
        "data/raw/study_table.parquet"
    output:
        "results/m2/meta_regression.txt"
    shell:
        "Rscript modules/m2_dag_meta/meta_regression.R"

rule m3_occupancy:
    input:
        "data/raw/detection_history.parquet"
    output:
        "results/m3/occupancy_summary.csv"
    shell:
        "Rscript modules/m3_occupancy/occupancy.R"

rule m1_mecor:
    input:
        "data/raw/study_table.parquet"
    output:
        "results/m1/mecor_report.txt"
    shell:
        "Rscript modules/m1_measurement_error/mecor_calibration.R"

rule m2_baggr:
    input:
        "data/raw/study_table.parquet"
    output:
        "results/m2/baggr_summary.txt"
    shell:
        "Rscript modules/m2_dag_meta/baggr_meta.R"

rule m3_shelf_priors:
    output:
        "results/m3/elicited_beta_priors.csv"
    shell:
        "Rscript modules/m3_occupancy/prior_elicitation_shelf.R"

rule m4_oncosimulr:
    output:
        "results/m4/oncosimulr_diversity.csv"
    shell:
        "Rscript modules/m4_dormancy/oncosimulr_selection.R"

rule m4_dormancy:
    input:
        "modules/m4_dormancy/model.py"
    output:
        "results/m4/dormancy_grid.csv"
    shell:
        "python modules/m4_dormancy/simulate.py"

rule m5_progression:
    input:
        "data/raw/first_progression.parquet"
    output:
        "results/m5/cif_brain_first_mut.csv",
        "results/m5/cif_brain_first_wt.csv"
    shell:
        "python modules/m5_progression/competing_risks.py"

rule m6_er_activity:
    input:
        "data/raw/expression_samples.parquet",
        "data/raw/expression.parquet",
        "data/raw/footprint_estrogen.csv"
    output:
        "results/m6/er_activity.csv"
    shell:
        "python modules/m6_er_activity/er_activity.py"
