# Copyright (c) 2019-2026 Arjan van der Velde, Weng Lab
# SPDX-License-Identifier: BSD-2-Clause

"""Parse and serialize ZDOCK and M-ZDOCK output files."""

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

FloatTriple = tuple[float, float, float]
IntTriple = tuple[int, int, int]


class ZDOCKError(ValueError):
    """Report invalid input or an operation unsupported by the parsed format."""


@dataclass
class Structure:
    """Describe a structure's filename, initial translation, and rotation."""

    filename: str = ""
    translation: FloatTriple = (0.0, 0.0, 0.0)
    rotation: FloatTriple = (0.0, 0.0, 0.0)

    def __repr__(self) -> str:
        """Return the structure in the metadata format used by ZDOCK."""
        return (
            f"{self.filename}\t{self.translation[0]:.3f}\t"
            f"{self.translation[1]:.3f}\t{self.translation[2]:.3f}"
        )


@dataclass
class Prediction:
    """Describe one predicted pose and its docking score."""

    translation: IntTriple = (0, 0, 0)
    rotation: FloatTriple = (0.0, 0.0, 0.0)
    score: float = 0.0
    ismzdock: bool = False

    def __repr__(self) -> str:
        """Return the prediction in its ZDOCK or M-ZDOCK line format."""
        if self.ismzdock:
            return (
                f"{self.rotation[0]:.6f}\t{self.rotation[1]:.6f}\t"
                f"{self.translation[0]}\t{self.translation[1]}\t{self.score:.2f}"
            )
        return (
            f"{self.rotation[0]:.6f}\t{self.rotation[1]:.6f}\t"
            f"{self.rotation[2]:.6f}\t{self.translation[0]}\t"
            f"{self.translation[1]}\t{self.translation[2]}\t{self.score:.3f}"
        )


class ZDOCK:
    """Parse one ZDOCK or M-ZDOCK output file.

    Args:
        filename: Path to the docking output file.
        n: Optional maximum number of predictions to retain. Zero retains none
            while still validating the file and parsing its metadata.
    """

    def __init__(self, filename: str | Path, n: int | None = None) -> None:
        if n is not None and n < 0:
            raise ValueError("Prediction limit cannot be negative")

        self._receptor = Structure()
        self._ligand = Structure()
        self._predictions: list[Prediction] = []
        self._boxsize = 0
        self._spacing = 0.0
        self._isswitched = False
        self._ismzdock = False
        self._isfixed = False
        self._version = 0
        self._symmetry = 0
        self._filename = Path(filename).expanduser()
        self._read(n)

    @property
    def receptor(self) -> Structure:
        """Return the receptor, or the sole structure for M-ZDOCK input."""
        return self._receptor

    @property
    def ligand(self) -> Structure:
        """Return the ligand for ZDOCK input.

        Raises:
            ZDOCKError: If the input is M-ZDOCK and has no separate ligand.
        """
        if self.ismzdock:
            raise ZDOCKError("ligand is not supported for M-ZDOCK")
        return self._ligand

    @property
    def npredictions(self) -> int:
        """Return the number of retained predictions."""
        return len(self._predictions)

    @property
    def predictions(self) -> tuple[Prediction, ...]:
        """Return retained predictions as an immutable sequence."""
        return tuple(self._predictions)

    @property
    def ismzdock(self) -> bool:
        """Return whether the input uses M-ZDOCK format."""
        return self._ismzdock

    @property
    def isfixed(self) -> bool:
        """Return whether the input uses the old fixed-receptor format."""
        return self._isfixed

    @property
    def symmetry(self) -> int:
        """Return the M-ZDOCK rotational symmetry.

        Raises:
            ZDOCKError: If the input is ordinary ZDOCK output.
        """
        if not self.ismzdock:
            raise ZDOCKError("symmetry is only supported for M-ZDOCK")
        return self._symmetry

    @property
    def boxsize(self) -> int:
        """Return the docking grid size."""
        return self._boxsize

    @property
    def filename(self) -> str:
        """Return the normalized input filename."""
        return str(self._filename)

    @property
    def version(self) -> int:
        """Return one for modern ZDOCK input and zero for legacy variants."""
        return self._version

    @property
    def spacing(self) -> float:
        """Return the docking grid spacing."""
        return self._spacing

    @property
    def lines(self) -> Iterator[str]:
        """Yield the normalized header and prediction lines."""
        return self._lines()

    @property
    def isswitched(self) -> bool:
        """Return whether receptor and ligand metadata were switched."""
        return self._isswitched

    @staticmethod
    def _process_line(line: str) -> list[str]:
        """Remove comments and split a nonempty line on arbitrary whitespace."""
        content = line.partition("#")[0].strip()
        return content.split() if content else []

    @staticmethod
    def _float_triple(fields: list[str]) -> FloatTriple:
        """Convert exactly three text fields to floating-point coordinates."""
        if len(fields) != 3:
            raise ZDOCKError("Expected three floating-point values")
        return (float(fields[0]), float(fields[1]), float(fields[2]))

    @staticmethod
    def _structure(fields: list[str]) -> Structure:
        """Convert filename and translation fields to a structure descriptor."""
        if len(fields) != 4:
            raise ZDOCKError("Expected a filename and three translation values")
        return Structure(
            filename=fields[0],
            translation=(float(fields[1]), float(fields[2]), float(fields[3])),
        )

    def _read(self, limit: int | None) -> None:
        """Read predictions and metadata from the configured file."""
        header: list[list[str]] = []
        header_done = False
        try:
            with self._filename.open("r", encoding="utf-8") as infile:
                for line_number, raw_line in enumerate(infile, start=1):
                    fields = self._process_line(raw_line)
                    if not fields:
                        continue
                    if len(fields) in (5, 7):
                        prediction = self._parse_prediction(fields, line_number)
                        self._predictions.append(prediction)
                        header_done = True
                        if limit is not None and 0 < limit <= len(self._predictions):
                            break
                    elif header_done:
                        raise ZDOCKError(f"Invalid prediction (line {line_number})")
                    else:
                        header.append(fields)
        except OSError as error:
            raise ZDOCKError(f"Unable to read {self._filename}: {error}") from error

        if limit == 0:
            self._predictions.clear()
        self._parse_header(header)

    def _parse_prediction(self, fields: list[str], line_number: int) -> Prediction:
        """Parse one prediction and enforce a consistent file format."""
        try:
            if len(fields) == 7:
                if self.ismzdock:
                    raise ZDOCKError(f"Invalid M-ZDOCK prediction (line {line_number})")
                return Prediction(
                    rotation=(float(fields[0]), float(fields[1]), float(fields[2])),
                    translation=(int(fields[3]), int(fields[4]), int(fields[5])),
                    score=float(fields[6]),
                )

            if not self.ismzdock and self._predictions:
                raise ZDOCKError(f"Invalid ZDOCK prediction (line {line_number})")
            self._ismzdock = True
            return Prediction(
                rotation=(float(fields[0]), float(fields[1]), 0.0),
                translation=(int(fields[2]), int(fields[3]), 0),
                score=float(fields[4]),
                ismzdock=True,
            )
        except ValueError as error:
            raise ZDOCKError(f"Invalid prediction values (line {line_number})") from error

    def _parse_header(self, header: list[list[str]]) -> None:
        """Parse format-specific metadata after predictions establish the format."""
        if self.ismzdock:
            self._parse_mzdock_header(header)
        else:
            self._parse_zdock_header(header)
        self._parse_structures(header)

    def _parse_zdock_header(self, header: list[list[str]]) -> None:
        """Parse modern or fixed-receptor ZDOCK grid metadata."""
        if len(header) == 5:
            self._version = 1
            self._isfixed = False
            expected_fields = 3
        elif len(header) == 4:
            self._version = 0
            self._isfixed = True
            expected_fields = 2
        else:
            raise ZDOCKError("ZDOCK header must have 4 or 5 rows")

        if len(header[0]) != expected_fields:
            raise ZDOCKError("Invalid ZDOCK grid metadata")
        try:
            self._boxsize = int(header[0][0])
            self._spacing = float(header[0][1])
            self._isswitched = bool(int(header[0][2])) if self._version else False
        except ValueError as error:
            raise ZDOCKError("Invalid ZDOCK grid metadata") from error

    def _parse_mzdock_header(self, header: list[list[str]]) -> None:
        """Parse M-ZDOCK grid and symmetry metadata."""
        if len(header) != 3 or len(header[0]) != 3:
            raise ZDOCKError("M-ZDOCK header must have 3 rows")
        try:
            self._boxsize = int(header[0][0])
            self._spacing = float(header[0][1])
            self._symmetry = int(header[0][2])
        except ValueError as error:
            raise ZDOCKError("Invalid M-ZDOCK grid metadata") from error
        if self._symmetry < 3:
            raise ZDOCKError("M-ZDOCK symmetry cannot be less than 3")

    def _parse_structures(self, header: list[list[str]]) -> None:
        """Parse initial rotations and translations for each structure."""
        try:
            if self.ismzdock:
                self._receptor = self._structure(header[2])
                self._receptor.rotation = self._float_triple(header[1])
                return

            if self.isfixed:
                self._receptor = self._structure(header[2])
                self._ligand = self._structure(header[3])
                self._ligand.rotation = self._float_triple(header[1])
                return

            receptor_index = 4 if self.isswitched else 3
            ligand_index = 3 if self.isswitched else 4
            self._receptor = self._structure(header[receptor_index])
            self._ligand = self._structure(header[ligand_index])
            self._receptor.rotation = self._float_triple(header[1])
            self._ligand.rotation = self._float_triple(header[2])
        except (IndexError, ValueError) as error:
            raise ZDOCKError("Invalid structure metadata") from error

    def _lines(self) -> Iterator[str]:
        """Yield normalized metadata followed by retained predictions."""
        if self.ismzdock:
            yield f"{self.boxsize}\t{self.spacing:.1f}\t{self.symmetry}"
            yield "\t".join(f"{value:.6f}" for value in self._receptor.rotation)
        elif self.isfixed:
            yield f"{self.boxsize}\t{self.spacing:.1f}"
            yield "\t".join(f"{value:.6f}" for value in self._ligand.rotation)
        else:
            yield f"{self.boxsize}\t{self.spacing:.1f}\t{int(self.isswitched)}"
            yield "\t".join(f"{value:.6f}" for value in self._receptor.rotation)
            yield "\t".join(f"{value:.6f}" for value in self._ligand.rotation)

        if self.ismzdock:
            yield repr(self._receptor)
        elif self.isswitched:
            yield repr(self._ligand)
            yield repr(self._receptor)
        else:
            yield repr(self._receptor)
            yield repr(self._ligand)
        yield from (repr(prediction) for prediction in self._predictions)

    def __repr__(self) -> str:
        """Return the complete normalized docking output."""
        return "\n".join(self._lines())
