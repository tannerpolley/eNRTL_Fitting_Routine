from __future__ import annotations

import json
import hashlib
import shutil
import tempfile
import unittest
from pathlib import Path

from packet_validator import PACKET_DIR, load_packet, validate_packet


class PacketValidator(unittest.TestCase):
    def copy_packet(self):
        tempdir = tempfile.TemporaryDirectory()
        packet = Path(tempdir.name) / "packet"
        shutil.copytree(PACKET_DIR, packet)
        self.addCleanup(tempdir.cleanup)
        return packet

    @staticmethod
    def refresh_payload_hash(packet):
        manifest_path = packet / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        payload_path = packet / "payload.json"
        for item in manifest["files"]:
            if item["path"] == "payload.json":
                content = payload_path.read_bytes()
                item["bytes"] = len(content)
                item["sha256"] = hashlib.sha256(content).hexdigest()
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    def test_positive_load_and_contract(self):
        manifest, payload = validate_packet()
        self.assertEqual(manifest["files"][0]["path"], "payload.json")
        self.assertEqual(
            payload["fit"]["observation_counts"],
            {"vle": 240, "speciation": 172, "calorimetry_fit": 66},
        )
        loaded = load_packet(
            expected_model="mea-enrtl-six-parameter-corrected-dh-v1",
            expected_species_order=manifest["model"]["species_order"],
            expected_reaction_order=manifest["model"]["reaction_order"],
            expected_property_basis="true",
            expected_domain="calibration_VLE",
            result_identity="baseline-six-parameter-physical-reaction-coordinates",
        )
        self.assertEqual(loaded["identity"], payload["identity"])

    def test_tampered_bytes_are_rejected(self):
        packet = self.copy_packet()
        payload = packet / "payload.json"
        payload.write_bytes(payload.read_bytes() + b"\n")
        with self.assertRaises(ValueError):
            validate_packet(packet)

    def test_semantic_tampering_is_rejected_after_hash_refresh(self):
        for case in (
            "model",
            "species",
            "reaction",
            "basis",
            "domain",
            "baseline",
            "unavailable",
            "counts",
        ):
            with self.subTest(case=case):
                packet = self.copy_packet()
                manifest_path = packet / "manifest.json"
                payload_path = packet / "payload.json"
                manifest = json.loads(manifest_path.read_text())
                payload = json.loads(payload_path.read_text())
                if case == "model":
                    manifest["model"]["identity"] = "wrong"
                elif case == "species":
                    manifest["model"]["species_order"] = ["wrong"]
                elif case == "reaction":
                    manifest["model"]["reaction_order"] = ["wrong"]
                elif case == "basis":
                    manifest["model"]["bases"]["property_basis"] = "apparent"
                elif case == "domain":
                    manifest["property_domains"]["calibration_VLE"]["loading"] = "wrong"
                elif case == "baseline":
                    payload["baseline_identity"] = "wrong"
                elif case == "unavailable":
                    manifest["unavailable_properties"].pop()
                else:
                    payload["fit"]["observation_counts"]["vle"] = 241
                manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
                payload_path.write_text(json.dumps(payload, indent=2) + "\n")
                self.refresh_payload_hash(packet)
                with self.assertRaises(ValueError):
                    validate_packet(packet)

    def test_non_finite_payload_is_rejected(self):
        packet = self.copy_packet()
        payload = packet / "payload.json"
        data = json.loads(payload.read_text())
        data["parameters"][0]["value"] = float("nan")
        payload.write_text(json.dumps(data))
        self.refresh_payload_hash(packet)
        with self.assertRaises(ValueError):
            validate_packet(packet)

    def test_declared_mismatches_are_rejected(self):
        cases = (
            {"expected_model": "wrong"},
            {"expected_species_order": ["wrong"]},
            {"expected_property_basis": "apparent"},
            {"expected_domain": "not-declared"},
            {"result_identity": "wrong"},
            {"requested_property": "total_reacting_solution_heat_capacity"},
        )
        for kwargs in cases:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                load_packet(**kwargs)


if __name__ == "__main__":
    unittest.main()
