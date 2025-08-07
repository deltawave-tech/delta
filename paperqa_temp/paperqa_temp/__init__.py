import warnings

# TODO: remove after refactoring all models to avoid using _* private vars
warnings.filterwarnings(
    "ignore", message="Valid config keys have changed in V2:", module="pydantic"
)


from paperqa_temp.agents import ask  # noqa: E402
from paperqa_temp.agents.main import agent_query  # noqa: E402
from paperqa_temp.agents.models import QueryRequest  # noqa: E402
from paperqa_temp.docs import Docs, PQASession, print_callback  # noqa: E402
from paperqa_temp.llms import (  # noqa: E402
    EmbeddingModel,
    HybridEmbeddingModel,
    LiteLLMEmbeddingModel,
    LiteLLMModel,
    LLMModel,
    LLMResult,
    NumpyVectorStore,
    SentenceTransformerEmbeddingModel,
    SparseEmbeddingModel,
    embedding_model_factory,
)
from paperqa_temp.settings import Settings, get_settings  # noqa: E402
from paperqa_temp.types import Answer, Context, Doc, DocDetails, Text  # noqa: E402
from paperqa_temp.version import __version__  # noqa: E402

__all__ = [
    "Answer",
    "Context",
    "Doc",
    "DocDetails",
    "Docs",
    "EmbeddingModel",
    "HybridEmbeddingModel",
    "LLMModel",
    "LLMResult",
    "LiteLLMEmbeddingModel",
    "LiteLLMModel",
    "NumpyVectorStore",
    "PQASession",
    "QueryRequest",
    "SentenceTransformerEmbeddingModel",
    "Settings",
    "SparseEmbeddingModel",
    "Text",
    "__version__",
    "agent_query",
    "ask",
    "embedding_model_factory",
    "get_settings",
    "print_callback",
]
