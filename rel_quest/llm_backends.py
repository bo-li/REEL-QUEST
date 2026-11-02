"""
LLM Backend abstraction for REL-QUEST AI Assistant.

Provides a common interface for different LLM providers:
  - OllamaBackend: local Ollama server
  - OpenAICompatibleBackend: OpenRouter, OpenAI, Groq, or any OpenAI-format API
"""

import json
import os
from typing import Any, Dict, Generator, List, Optional

import requests

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OLLAMA_BASE_URL = os.environ.get("OLLAMA_BASE_URL", "http://localhost:11434")
REQUEST_TIMEOUT = 120

OPENROUTER_BASE_URL = "https://openrouter.ai/api/v1"


def get_config(name: str, default: str = "") -> str:
    """
    Read a setting from Streamlit secrets, falling back to the environment.

    Streamlit Community Cloud supplies secrets through ``st.secrets``; local
    runs and Docker use environment variables. Accessing ``st.secrets`` raises
    when no secrets file exists, so the whole lookup is guarded and the
    streamlit import is deferred — this module must stay importable from
    plain scripts, the REST API and pytest.
    """
    try:
        import streamlit as st

        value = st.secrets.get(name)
        if value:
            return str(value)
    except Exception:
        pass
    return os.environ.get(name, default)


# Zero-cost OpenRouter IDs that carry no ``:free`` suffix. ``openrouter/free``
# is the "Free Models Router", which dispatches to whichever free models are
# currently up — priced at 0 for both prompt and completion.
FREE_ROUTER_IDS = frozenset({"openrouter/free"})


def is_free_model(model: str) -> bool:
    """
    True for OpenRouter model IDs that cost nothing to call.

    That is the ``:free`` suffix, plus the suffix-less free routers in
    ``FREE_ROUTER_IDS``. The shared owner-supplied key must never be pointed at
    a paid model, so every call site that uses it checks this first.
    """
    if not model:
        return False
    model = model.strip()
    return model.endswith(":free") or model in FREE_ROUTER_IDS


# ---------------------------------------------------------------------------
# Base class
# ---------------------------------------------------------------------------

class LLMBackend:
    """Abstract base class for LLM backends."""

    def chat_stream(self, messages: List[Dict]) -> Generator[str, None, None]:
        """Stream response tokens as strings."""
        raise NotImplementedError

    def chat_sync(self, messages: List[Dict], tools: Optional[List] = None) -> Dict:
        """Synchronous chat call. Returns a response dict with 'message' key."""
        raise NotImplementedError

    def is_available(self) -> bool:
        """Return True if the backend can be reached."""
        raise NotImplementedError


# ---------------------------------------------------------------------------
# Ollama backend
# ---------------------------------------------------------------------------

class OllamaBackend(LLMBackend):
    """Local Ollama server backend."""

    def __init__(self, base_url: str = OLLAMA_BASE_URL, model: str = "qwen3.5:0.8b"):
        self.base_url = base_url.rstrip("/")
        self.model = model

    def chat_sync(self, messages: List[Dict], tools: Optional[List] = None) -> Dict:
        payload: Dict[str, Any] = {
            "model": self.model,
            "messages": messages,
            "stream": False,
            "options": {"temperature": 0.1},
        }
        if tools:
            payload["tools"] = tools
        resp = requests.post(
            f"{self.base_url}/api/chat",
            json=payload,
            timeout=REQUEST_TIMEOUT,
        )
        resp.raise_for_status()
        return resp.json()

    def chat_stream(self, messages: List[Dict]) -> Generator[str, None, None]:
        payload: Dict[str, Any] = {
            "model": self.model,
            "messages": messages,
            "stream": True,
            "options": {"temperature": 0.1},
        }
        with requests.post(
            f"{self.base_url}/api/chat",
            json=payload,
            stream=True,
            timeout=REQUEST_TIMEOUT,
        ) as resp:
            resp.raise_for_status()
            for line in resp.iter_lines():
                if line:
                    try:
                        chunk = json.loads(line)
                        content = chunk.get("message", {}).get("content", "")
                        if content:
                            yield content
                        if chunk.get("done", False):
                            break
                    except json.JSONDecodeError:
                        continue

    def is_available(self) -> bool:
        try:
            r = requests.get(f"{self.base_url}/api/tags", timeout=5)
            return r.status_code == 200
        except Exception:
            return False

    def list_models(self) -> List[str]:
        try:
            r = requests.get(f"{self.base_url}/api/tags", timeout=5)
            r.raise_for_status()
            return [m["name"] for m in r.json().get("models", [])]
        except Exception:
            return []

    def pull_model(self, name: str) -> bool:
        try:
            resp = requests.post(
                f"{self.base_url}/api/pull",
                json={"name": name, "stream": False},
                timeout=REQUEST_TIMEOUT,
            )
            return resp.status_code == 200
        except Exception:
            return False


# ---------------------------------------------------------------------------
# OpenAI-compatible backend (OpenRouter, OpenAI, Groq, etc.)
# ---------------------------------------------------------------------------

class OpenAICompatibleBackend(LLMBackend):
    """
    Backend for any OpenAI /v1/chat/completions compatible API.

    Works with:
      - OpenRouter (https://openrouter.ai/api/v1)
      - OpenAI    (https://api.openai.com/v1)
      - Groq      (https://api.groq.com/openai/v1)
      - Local vLLM / LiteLLM proxy

    Requires the `openai` Python package (pip install openai>=1.0.0).
    """

    def __init__(
        self,
        api_key: str,
        model: str = "openrouter/free",
        base_url: str = OPENROUTER_BASE_URL,
    ):
        try:
            import openai as _openai
        except ImportError as exc:
            raise ImportError(
                "The `openai` package is required for cloud LLM backends. "
                "Install it with: pip install openai>=1.0.0"
            ) from exc

        self.model = model
        self._client = _openai.OpenAI(api_key=api_key, base_url=base_url)

    def chat_sync(self, messages: List[Dict], tools: Optional[List] = None) -> Dict:
        """Returns an Ollama-style dict so agent.py call sites are unchanged."""
        resp = self._client.chat.completions.create(
            model=self.model,
            messages=messages,
            stream=False,
            temperature=0.1,
        )
        content = resp.choices[0].message.content or ""
        return {"message": {"role": "assistant", "content": content}}

    def chat_stream(self, messages: List[Dict]) -> Generator[str, None, None]:
        stream = self._client.chat.completions.create(
            model=self.model,
            messages=messages,
            stream=True,
            temperature=0.1,
        )
        for chunk in stream:
            delta = chunk.choices[0].delta.content
            if delta:
                yield delta

    def is_available(self) -> bool:
        try:
            self._client.models.list()
            return True
        except Exception:
            return False
