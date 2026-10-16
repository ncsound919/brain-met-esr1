"""
Pydantic schemas for the shared data spine.
Every module reads from these tables.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, model_validator

Compartment = Literal["brain", "liver", "bone", "lung", "lymph_node", "other"]
DetectionMethod = Literal["tumor_tissue_NGS", "ctDNA_NGS", "ctDNA_ddPCR", "WES", "IHC", "other"]
DenominatorUnit = Literal["patients", "samples"]
Esr1Status = Literal["mut", "wt"]
ProgressionEvent = Literal["brain_first", "nonCNS_first", "censored"]


class DetectionHistory(BaseModel):
    """
    One row per patient, compartment, timepoint, and method.
    'detected=False' means nondetection, not absence.
    """

    patient_id: str
    compartment: Compartment
    timepoint_days: int = Field(ge=0)
    detection_method: DetectionMethod
    assay_platform: str
    limit_of_detection: float | None = Field(default=None, ge=0.0)
    detected: bool
    ctdna_fraction: float | None = Field(default=None, ge=0.0, le=1.0)
    denominator_unit: DenominatorUnit

    @model_validator(mode="after")
    def check_ctdna_fraction(self) -> DetectionHistory:
        if "ctDNA" in self.detection_method and self.ctdna_fraction is None:
            raise ValueError(
                "ctDNA detection rows must include ctdna_fraction; "
                "otherwise nondetection is ambiguous."
            )
        return self


class StudyRecord(BaseModel):
    """One row per included study."""

    study_id: str
    citation: str
    subtype_scope: str
    n_brain: int = Field(ge=0)
    n_nonCNS: int = Field(ge=0)
    esr1_brain_mut: int = Field(ge=0)
    esr1_nonCNS_mut: int = Field(ge=0)
    detection_method: DetectionMethod
    assay_platform: str
    denominator_unit: DenominatorUnit
    reported_or: float | None = None
    reported_or_lo: float | None = None
    reported_or_hi: float | None = None

    @model_validator(mode="after")
    def check_counts(self) -> StudyRecord:
        if self.esr1_brain_mut > self.n_brain:
            raise ValueError(f"{self.study_id}: esr1_brain_mut exceeds n_brain")
        if self.esr1_nonCNS_mut > self.n_nonCNS:
            raise ValueError(f"{self.study_id}: esr1_nonCNS_mut exceeds n_nonCNS")
        return self


class FirstProgression(BaseModel):
    """One row per patient: first progression site (competing event) or censoring.

    M5 treats brain_first vs nonCNS_first as competing risks, not censoring:
    standard KM treating the other site as censored biases the absolute risk.
    """

    patient_id: str
    esr1_status: Esr1Status
    time_days: int = Field(ge=0)
    event: ProgressionEvent
    detection_method: DetectionMethod
    subtype_scope: str = "ER+/HER2-"


class ExpressionSample(BaseModel):
    """One row per RNA sample. The matrix lives in expression.parquet
    (tidy: sample_id, gene, log_norm_count). M6 joins on sample_id."""

    sample_id: str
    patient_id: str
    compartment: Compartment
    esr1_status: Esr1Status
    assay_platform: str


class ExpressionRow(BaseModel):
    """One tidy expression row: sample_id, gene (HGNC), log_norm_count."""

    sample_id: str
    gene: str
    log_norm_count: float


def check_duplicates(det_rows: list[dict]) -> None:
    """Fail on contradictory duplicate detection rows.

    The same (patient_id, compartment, timepoint_days, detection_method)
    must not appear twice with different `detected` values.
    """
    seen: dict[tuple, bool] = {}
    for row in det_rows:
        key = (
            row["patient_id"],
            row["compartment"],
            row["timepoint_days"],
            row["detection_method"],
        )
        if key in seen and seen[key] != row["detected"]:
            raise ValueError(
                f"Contradictory duplicate detection row for {key}: {seen[key]} vs {row['detected']}"
            )
        seen[key] = row["detected"]
