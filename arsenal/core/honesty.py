#!/usr/bin/env python3
"""
honesty.py -- the trust layer for the RH arsenal.

Three tiers, never conflated (Charisma Report/Handle + the project's proof-pointer rule):

    CLAIMED  -- asserted, no evidence. Default for anything unrun.
    REPORT   -- describable/analyzable now: a scorer ran, a structural check passed.
                Real work, but NOT proof.
    HANDLE   -- earned ONLY by a real run producing a real artifact that exists on disk.
                assert_handle() REFUSES to stamp HANDLE if the proof-pointer file is absent.

This module makes false certification structurally hard: you cannot construct a Result at
HANDLE tier without pointing at a file that exists, and the banned-phrase guard blocks any
tool from emitting "proven RH" / "certified RH" style overclaims.

No RH tool in this arsenal is allowed to claim it proved or certified the Riemann Hypothesis.
Finite/numeric results must say "finite" / "numeric" in the same object.
"""
from __future__ import annotations
import os, json, hashlib
from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Any, Optional


class Tier(str, Enum):
    CLAIMED = "CLAIMED"   # asserted, unverified
    REPORT = "REPORT"     # analyzed/ran, not proof
    HANDLE = "HANDLE"     # earned by a real artifact on disk

    def __str__(self):  # so it prints as the bare word
        return self.value


# Phrases a tool must never emit as a claim about RH itself.
_BANNED = [
    "proven rh", "proved rh", "proves rh", "proving rh", "rh is proven",
    "rh proven", "rh is proved", "certified rh", "rh certified",
    "riemann hypothesis is proven", "riemann hypothesis proved",
    "proves the riemann", "proves riemann", "solves rh", "rh solved",
    "solved rh", "solves the riemann", "proof of rh", "proof of the riemann",
    "complete proof of rh", "proof of rh complete", "qed rh",
    "rh is true (proven)", "we have proven the riemann", "proof of riemann hypothesis",
]


class OverclaimError(Exception):
    pass


class NoEvidenceError(Exception):
    pass


def guard_text(text: str) -> str:
    """Raise if a claim string overclaims a proof of RH. Returns text unchanged if clean."""
    low = text.lower()
    for b in _BANNED:
        if b in low:
            raise OverclaimError(
                f"Refusing to emit banned RH-proof overclaim: {b!r} in {text!r}. "
                f"Finite/numeric results are not a proof of RH."
            )
    return text


@dataclass
class Result:
    """A single honest result. Tier is enforced against evidence at construction time."""
    claim: str                     # what is asserted
    tier: Tier                     # CLAIMED / REPORT / HANDLE
    method: str                    # how it was obtained
    proof_pointer: Optional[str] = None   # path to an artifact that backs a HANDLE claim
    data: dict = field(default_factory=dict)
    tool: str = ""
    finite: bool = True            # is this a finite/numeric (non-RH-closing) result?

    def __post_init__(self):
        if isinstance(self.tier, str):
            self.tier = Tier(self.tier)
        guard_text(self.claim)
        if self.tier == Tier.HANDLE:
            if not self.proof_pointer:
                raise NoEvidenceError(
                    f"HANDLE tier requires a proof_pointer. claim={self.claim!r}")
            if not os.path.exists(self.proof_pointer):
                raise NoEvidenceError(
                    f"HANDLE tier proof_pointer does not exist on disk: "
                    f"{self.proof_pointer!r} for claim={self.claim!r}")

    def to_dict(self) -> dict:
        d = asdict(self)
        d["tier"] = str(self.tier)
        return d

    def line(self) -> str:
        p = f"  [proof: {self.proof_pointer}]" if self.proof_pointer else ""
        fin = "" if not self.finite else "  (finite/numeric)"
        return f"[{self.tier}] {self.claim}{fin}{p}"


def assert_handle(condition: bool, claim: str, proof_pointer: str, method: str,
                  tool: str = "", data: dict | None = None) -> Result:
    """
    Construct a HANDLE result ONLY if `condition` is true AND the proof_pointer file exists.
    Anything else downgrades to REPORT (condition false) or raises (no evidence file).
    This is the single choke point through which "it works" claims must pass.
    """
    if not condition:
        return Result(claim="NOT ESTABLISHED: " + claim, tier=Tier.REPORT,
                      method=method, tool=tool, data=data or {})
    return Result(claim=claim, tier=Tier.HANDLE, method=method,
                  proof_pointer=proof_pointer, tool=tool, data=data or {})


def sha256_of(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def digest_obj(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]


def selftest(verbose=True) -> bool:
    ok = True
    # 1. banned phrase is blocked
    try:
        Result(claim="this proves RH", tier=Tier.REPORT, method="x")
        ok = False; print("FAIL: overclaim not blocked")
    except OverclaimError:
        pass
    # 2. HANDLE without evidence is blocked
    try:
        Result(claim="finite margin positive", tier=Tier.HANDLE, method="x",
               proof_pointer="/no/such/file.json")
        ok = False; print("FAIL: HANDLE without evidence not blocked")
    except NoEvidenceError:
        pass
    # 3. false condition downgrades to REPORT
    r = assert_handle(False, "matrix is PD", "/tmp/x", "ldlt")
    if r.tier != Tier.REPORT:
        ok = False; print("FAIL: false condition did not downgrade")
    # 4. real evidence file yields HANDLE
    tmp = "/tmp/_honesty_selftest_evidence.json"
    with open(tmp, "w") as f:
        json.dump({"pivots_all_positive": True}, f)
    r = assert_handle(True, "finite 2x2 matrix PD over Q", tmp, "exact-ldlt")
    if r.tier != Tier.HANDLE:
        ok = False; print("FAIL: real evidence did not yield HANDLE")
    os.remove(tmp)
    print("honesty.selftest:", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    selftest()
