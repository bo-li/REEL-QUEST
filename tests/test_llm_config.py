"""
Tests for the LLM backend configuration helpers.

These guard the two pieces the Streamlit Cloud deployment depends on:
`get_config`, which must survive a missing secrets file and fall back to the
environment, and `is_free_model`, the predicate that keeps the app owner's
shared OpenRouter key from ever being pointed at a paid model.
"""

import os
import unittest
from unittest import mock

from rel_quest.llm_backends import get_config, is_free_model


class TestIsFreeModel(unittest.TestCase):

    def test_accepts_free_tier_ids(self):
        for model in ("nvidia/nemotron-3-super-120b-a12b:free",
                      "nvidia/nemotron-3-nano-omni-30b-a3b-reasoning:free",
                      "  minimax/minimax-m3:free  "):
            self.assertTrue(is_free_model(model), model)

    def test_accepts_suffixless_free_routers(self):
        # openrouter/free is priced at zero but carries no ":free" suffix.
        self.assertTrue(is_free_model("openrouter/free"))
        self.assertTrue(is_free_model("  openrouter/free  "))

    def test_rejects_paid_and_empty(self):
        for model in ("openai/gpt-4o-mini", "anthropic/claude-opus-5",
                      "openrouter/auto", "minimax/minimax-m3", "qwen2.5:7b",
                      "", None):
            self.assertFalse(is_free_model(model), model)

    def test_free_suffix_must_end_the_id(self):
        # ":free" elsewhere in the string is not a free-tier model.
        self.assertFalse(is_free_model("vendor/model:free-preview"))


class TestGetConfig(unittest.TestCase):
    """
    No secrets.toml exists during a plain pytest run, so `st.secrets` raises —
    exercising the fallback path that Cloud-free deployments rely on.
    """

    def test_reads_environment(self):
        with mock.patch.dict(os.environ, {"RELQUEST_TEST_KEY": "from-env"}):
            self.assertEqual(get_config("RELQUEST_TEST_KEY"), "from-env")

    def test_returns_default_when_unset(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertEqual(get_config("RELQUEST_ABSENT", "fallback"), "fallback")
            self.assertEqual(get_config("RELQUEST_ABSENT"), "")

    def test_secrets_take_precedence_over_environment(self):
        fake_st = mock.MagicMock()
        fake_st.secrets.get.return_value = "from-secrets"
        with mock.patch.dict(os.environ, {"RELQUEST_TEST_KEY": "from-env"}), \
                mock.patch.dict("sys.modules", {"streamlit": fake_st}):
            self.assertEqual(get_config("RELQUEST_TEST_KEY"), "from-secrets")

    def test_falls_back_when_secrets_raise(self):
        fake_st = mock.MagicMock()
        fake_st.secrets.get.side_effect = FileNotFoundError("no secrets.toml")
        with mock.patch.dict(os.environ, {"RELQUEST_TEST_KEY": "from-env"}), \
                mock.patch.dict("sys.modules", {"streamlit": fake_st}):
            self.assertEqual(get_config("RELQUEST_TEST_KEY"), "from-env")


if __name__ == "__main__":
    unittest.main()
