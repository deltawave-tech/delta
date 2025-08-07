from dotenv import load_dotenv


load_dotenv()

from enum import Enum
from typing import Dict, Union, List, Optional, Any, Tuple
from pydantic import BaseModel
from anthropic import Anthropic
import anthropic
import copy
import re
import base64
import logging
import time
import random
from src.tools.tool_definitions import tool_registry
## new agent file

# Add litellm import
import litellm
import json
class ClaudeResponse(BaseModel):
    content: str
    content_blocks: Optional[List[Dict[str, Any]]] = None
    history: Optional[List[Dict[str, Any]]] = None
    results: Optional[List[Any]] = None
    function_call: Optional[Any] = None
    price: Optional[float] = None
    class Config: arbitrary_types_allowed = True

class GeminiResponse(BaseModel): # Not directly used by the new Agent but kept for context
    content: str
    history: Optional[List[Dict[str, Any]]] = None
    results: Optional[List[Any]] = None
    function_call: Optional[Any] = None
    price: Optional[float] = None
provider2model = {"sonnet-4":"claude-sonnet-4-20250514",
                  "sonnet-3.7":"claude-3-7-sonnet-20250219",
                  "o3":"o3-2025-04-16",
                  "gpt4.1": "gpt-4.1-2025-04-14",
                  "gemini":"gemini/gemini-2.5-pro-preview-06-05",
                  "qwen": "openrouter/qwen/qwen3-235b-a22b",
                  "deepseek":"deepseek/deepseek-chat-v3-0324",
                  "mistral": "openrouter/mistralai/magistral-medium-2506"}


def transform_to_openai_format(input_dict):
    """
    Transforms a dictionary from a custom tool format to the OpenAI function call definition format.

    Args:
        input_dict (dict): A dictionary with 'name', 'description', and 'input_schema'.
                           'input_schema' should itself be a dictionary containing
                           'type', 'properties', and optionally 'required'.

    Returns:
        dict: A dictionary in the OpenAI function call definition format.
    """
    if not isinstance(input_dict, dict):
        raise TypeError("Input must be a dictionary.")

    if not all(key in input_dict for key in ["name", "description", "input_schema"]):
        raise ValueError("Input dictionary must contain 'name', 'description', and 'input_schema' keys.")

    if not isinstance(input_dict.get("input_schema"), dict):
        raise ValueError("'input_schema' must be a dictionary.")

    openai_function_definition = {
        "type": "function",
        "function": {
            "name": input_dict["name"],
            "description": input_dict["description"],
            "parameters": input_dict["input_schema"]  # Directly assign input_schema
        }
    }
    return openai_function_definition

class Agent:
    def __init__(
        self,
        title: str,
        expertise: str,
        goal: str,
        role:str,
        model: str = None,
        tools: Optional[List[dict]] = None,
        max_tokens: int = 4096*3,
        temperature: float = 0.1,
        llm_provider = None
    ) -> None:
        self.title = title
        self.expertise = expertise
        self.role = role
        self.goal = goal
        self.model = model
        self.max_tokens = max_tokens
        self.temperature = temperature
        self.tools = tools or []
        self.tool_history = [] # For storing history_with_tools at the end of successful conversations
        if model is None: 
            self.model = provider2model[llm_provider]


# "parent_molecule_id": "ID of the parent molecule if it exists, Database Agent and AI Expert will not have a parent molecule id",
    @property
    def prompt(self) -> str:
        molecule_format_guide = """
        When presenting molecules, always use the following structured format:
        <molecules>
        {
            "current_iteration": 1,
            "molecules": [
                {
                    "smiles": "SMILES_STRING",
                    "friendly_id": "friendly id of molecule eg: DA:I1:N0:G0",
                    "source": "(parent molecule id of modified molecule)"
                    "cycle": 1,
                    "status": "(active, inactive, de_novo, modified)",
                    "iteration": iteration number,
                    "agent": "AGENT_NAME",
                    "metrics": {
                        "docking": "DOCKING_SCORE",
                        "qed": "QED_SCORE",
                        "sa_score": "SAS_SCORE",
                        "logp": "LOG_P",
                        "molecular_weight": "MOLECULAR_WEIGHT",
                        "plip_interactions": "(plip interactions for the molecule)"
                    },
                    "rationale": "(rationale for the molecule)"
                }
            ],
            "summaries": [
                {
                    "iteration": iteration number1,
                    "summary": "SUMMARY"
                }
            ]
        }
        </molecules>
                """
        molecule_format_guide_simplified = """
        "molecules": [
                {
                    "smiles": "SMILES_STRING",
                    "friendly_id": "friendly id of molecule eg: DA:I1:N0:G0",
                    "source": "(parent molecule id of modified molecule)"
                    "agent": "AGENT_NAME",
                    "rationale": "(rationale for the molecule)"
                }
            ],
        }
        </molecules>
        """
        molecule_id_format_guide = """
Druing the discovery process, molecules will be presented with a friendly_id. This friendly_id is a unique identifier that tracks its complete history and lineage through our optimization process to help to follow the SAR and parent-child relationships.

## Molecule ID System - Tracking Molecular Evolution

Each molecule receives a unique identifier that tracks its complete history and lineage through our optimization process.

**Format:** `AgentCode:I<iteration>:N<position>:G<generation>`

This ID (friendly_id) system serves as a molecular genealogy, allowing us to:
- Track which agent created or modified each molecule
- Understand the evolutionary path from original to optimized compounds
- Maintain parent-child relationships across modifications

**Components:**
- **AgentCode**: Two-letter code representing the creator agent
  - DA = Database Agent
  - AI = AI Expert
  - MC = Medicinal Chemist
  - LA = Literature Agent
  - RA = Ranking Agent
  - PR = Principal Researcher
  - SC = Scientific Critic
  - SP = Summary Parser

- **I<iteration>**: The research iteration number when molecule was introduced (e.g., I1, I2, I3)
- **N<position>**: The molecule's position in the agent's output list, starting from N1 (resets per agent output)
- **G<generation>**: The modification depth showing molecular evolution
  - G0 = original molecule (from database or de novo)
  - G1 = first modification
  - G2 = second modification
  - And so on...

**Example Evolution:**
1. Database retrieves: `DA:I1:N3:G0` (original from ChEMBL, third molecule in list, iteration 1)
2. Medicinal Chemist modifies it: `MC:I2:N1:G1` (first modification, first molecule in chemist's list, iteration 2)
3. Further optimization: `MC:I3:N2:G2` (second modification, second molecule in list, iteration 3)

This allows us to trace any molecule back to its origins and understand its optimization journey.
If a molecule is a child of a parent molecule, the parent molecule's friendly_id must be included in the molecule's source field
"""
        prompt_medchem = '''## Special Instructions - Molecule Modification Protocol for you medicinal chemist agent

Your primary role is to **modify existing molecules** to improve their properties. You MUST always specify the parent molecule(s) that inspired your modification to maintain the evolutionary chain and track Structure-Activity Relationships (SAR).

**Modification Rules & Parentage:**

1.  **ALWAYS Identify Parent(s):**
    *   When proposing a new molecule, you **MUST** state which molecule(s) you are modifying or were inspired by. Provide the full `friendly_id` of the parent(s).
    *   **Single Parent:** If one primary molecule was modified, provide its `friendly_id` as the `source`.
        *   Example: `source: "AI:I1:N4:G0"`
    *   **Multiple Parents/Inspirations:** If your new molecule is a hybrid, a significant departure, or inspired by multiple existing molecules, provide a list of `friendly_id`s as the `source`.
        *   Example: `source: ["AI:I1:N2:G0", "DA:I1:N5:G0"]`
    *   **"Most Inspiring" Parent:** Even if your modification is substantial, choose the `friendly_id` of the molecule that most significantly influenced your design. **Do NOT label your own outputs as "de novo".** All your work should build upon or be inspired by prior molecules in the history.

2.  **Generation Tracking (Automatic):**
    *   The system automatically assigns the new generation number (G-number) based on the parent(s).
        *   Parent G0 → Child becomes G1
        *   Parent G1 → Child becomes G2, etc.
    Important note: a child molecule can 

3.  **Output Format for Modified Molecules:**
    *   Your output for each molecule should follow the established JSON structure, paying close attention to the `source` field.
    *   **Example (Single Parent):**
        ```json
        {
            "smiles": "O=C(C1=CC(Br)=CNC1)n1ccc2ccc3ncccc3c21",
            "friendly_id": "MC:I1:N4:G1", // Automatically assigned based on your iteration & agent
            "source": "AI:I1:N4:G0",     // YOU MUST PROVIDE THIS
            "status": "modified",
            // ... other fields like metrics, rationale ...
        }
        ```
    *   **Example (Multiple Parents/Inspirations):**
        ```json
        {
            "smiles": "NEW_HYBRID_SMILES",
            "friendly_id": "MC:I1:N5:G1",
            "source": ["AI:I1:N2:G0", "DA:I1:N5:G0"], // List of inspiring parent IDs
            "status": "modified_hybrid", // Or just "modified"
            // ... other fields ...
        }
        ```

**Example Workflow:**

*   You decide to modify a molecule with `friendly_id`: `AI:I1:N4:G0` (Generation G0).
*   Your modification is made in Iteration 1 (I1) by you, the Medicinal Chemist (MC).
*   Your new molecule will be assigned a `friendly_id` like `MC:I1:N_new:G1`.
*   You **MUST** set `source: "AI:I1:N4:G0"` in your output.

*   If you then take `MC:I1:N_new:G1` and further modify it in Iteration 2:
*   The next molecule will be `MC:I2:N_another:G2`.
*   You **MUST** set `source: "MC:I1:N_new:G1"` for this second modification.

**CRITICAL IMPORTANCE:**
*   **No "De Novo" from You:** As the Medicinal Chemist, your role is to optimize and evolve. Avoid using "De novo design..." or similar phrases in the `source` field. Trace every molecule back to one or more parents from the existing pool.
*   **Complete Parent Information:** Without specifying the parent `friendly_id`(s) in the `source` field, the system cannot track molecular evolution. This breaks the SAR analysis and obscures the optimization pathway.

Focus on iterating upon existing structures, drawing inspiration from their features, and clearly documenting the lineage of your new designs through the `source` field.
! Keep in mind, when a source molecule has G={ i }, the child molecule will necessarily have G={i+1}!
If a molecule is a child of a parent molecule, the parent molecule's friendly_id must be included in the molecule's source field
'''
        include_molecule_guide = self.title in ["Database Agent", "AI Expert", "Medicinal Chemist"]
        molecule_guide_section = molecule_format_guide if include_molecule_guide else molecule_format_guide_simplified
        if self.title == 'Summary Parser':
            molecule_guide_section = ""
        prompt_medchem = prompt_medchem if self.title == "Medicinal Chemist" else ""
        return (
            f"You are an AI assistant acting as {self.title}. "
            f"Your expertise is in {self.expertise}. "
            f"Your goal is to {self.goal}. "
            f"Your role is to {self.role}. "
            f"Please maintain this role throughout our conversation.\n\n"
            f"{molecule_guide_section}"
            f"{molecule_id_format_guide}"
            f"{prompt_medchem}"
        )
    def __hash__(self) -> int:
        return hash(self.title)

    def __eq__(self, other: object) -> bool:
        
        if not isinstance(other, Agent): return False
        return (
            self.title == other.title and self.expertise == other.expertise and
            self.goal == other.goal and self.role == other.role and
            self.model == other.model # Compare model string
        )

    def __str__(self) -> str:
        return self.title

    def _get_completion(self, history: List[Dict[str, Any]],
                        temperature: float,
                        max_tokens: int,
                        logger: logging.Logger,
                        max_retries: int = 5,
                        base_delay: float = 1.0) -> Tuple[litellm.ModelResponse, float]:
        messages_with_system_prompt = []
        if self.prompt:
            messages_with_system_prompt.append({"role": "system", "content": self.prompt})
        messages_with_system_prompt.extend(history)

        # O3 models only support temperature=1
        if "o3" in self.model.lower():
            temperature = 1.0

        kwargs = {
            "model":self.model,
            # "model":'gpt-4o-mini-2024-07-18',
            # "model":'gpt-4.1-2025-04-14',
            # "model":'claude-opus-4-20250514',
            # "model":'claude-sonnet-4-20250514',
            # "model":'gemini/gemini-2.5-pro-preview-06-05',
            "messages": messages_with_system_prompt,
            "temperature": temperature,
            "max_tokens": max_tokens,
        }
        
        # Add Gemini-specific safety settings to prevent content blocking
        if "gemini" in self.model.lower():
            kwargs["safety_settings"] = [
                {
                    "category": "HARM_CATEGORY_HARASSMENT",
                    "threshold": "BLOCK_NONE"
                },
                {
                    "category": "HARM_CATEGORY_HATE_SPEECH", 
                    "threshold": "BLOCK_NONE"
                },
                {
                    "category": "HARM_CATEGORY_SEXUALLY_EXPLICIT",
                    "threshold": "BLOCK_NONE"
                },
                {
                    "category": "HARM_CATEGORY_DANGEROUS_CONTENT",
                    "threshold": "BLOCK_NONE"
                }
            ]
        
        if self.tools:
            kwargs["tools"] = [transform_to_openai_format(tool) for tool in self.tools]
            logger.info(f"DEBUG: tools: {kwargs['tools']}")
            kwargs["tool_choice"] = "auto"
        # logger.info(f"DEBUG: Sending to litellm.completion with kwargs: {json.dumps(kwargs, indent=2)}")
        logger.info(f'chat history for chat completion: {json.dumps(messages_with_system_prompt, indent=2)}')
        # logger.info(f'system prompt: {self.prompt}')
        
        last_exception = None
        
        for attempt in range(max_retries + 1):
            try:
                response: litellm.ModelResponse = litellm.completion(**kwargs)
                try:
                    price = litellm.completion_cost(completion_response=response)
                except Exception as cost_error:
                    logger.warning(f"Could not calculate completion cost for model {self.model}: {str(cost_error)}")
                    price = 0.0
                
                # If we succeed after retries, log the success
                if attempt > 0:
                    logger.info(f"Successfully completed request after {attempt} retries")
                
                return response, price
                
            except (litellm.exceptions.InternalServerError, 
                    litellm.exceptions.RateLimitError,
                    litellm.exceptions.ServiceUnavailableError,
                    litellm.exceptions.APIError) as e:
                
                last_exception = e
                
                # Check if it's an Anthropic overload error or other retryable error
                is_retryable = (
                    "overloaded_error" in str(e).lower() or 
                    "overloaded" in str(e).lower() or
                    "rate limit" in str(e).lower() or
                    "service unavailable" in str(e).lower() or
                    "timeout" in str(e).lower() or
                    "429" in str(e) or  # Rate limit HTTP status
                    "529" in str(e) or  # Service overloaded HTTP status
                    "503" in str(e)     # Service unavailable HTTP status
                )
                
                if not is_retryable or attempt == max_retries:
                    if attempt == max_retries:
                        logger.error(f"Max retries ({max_retries}) exceeded for completion request")
                    else:
                        logger.error(f"Non-retryable error during completion: {str(e)}")
                    raise e
                
                # Calculate delay with exponential backoff and jitter
                delay = base_delay * (2 ** attempt) + random.uniform(0, 1)
                logger.warning(f"Retryable error on attempt {attempt + 1}/{max_retries + 1}: {str(e)}")
                logger.info(f"Retrying in {delay:.2f} seconds...")
                time.sleep(delay)
                
            except (litellm.exceptions.BudgetExceededError,
                    litellm.exceptions.InvalidRequestError,
                    litellm.ContextWindowExceededError) as e:
                # These are non-retryable errors
                logger.error(f"Non-retryable error during completion: {str(e)}")
                raise
                
            except Exception as e:
                last_exception = e
                
                # For unknown exceptions, be more conservative and only retry a few times
                if attempt < min(2, max_retries):
                    delay = base_delay * (2 ** attempt) + random.uniform(0, 1)
                    logger.warning(f"Unknown error on attempt {attempt + 1}, retrying in {delay:.2f} seconds: {str(e)}")
                    time.sleep(delay)
                else:
                    logger.error(f"Unexpected error during completion after {attempt} retries: {str(e)}")
                    raise
        
        # This should never be reached, but just in case
        if last_exception:
            raise last_exception
        else:
            raise Exception("Unknown error: maximum retries exceeded without specific exception")

    def format_message(self, role: str, content: Union[str, List[Dict[str, Any]], Dict[str, Any]]):
        # This format is for Claude-style content blocks, used for the `history` in `ClaudeResponse`
        # and for initial user message. LiteLLM messages can also be simpler.
        if isinstance(content, str):
            return {"role": role, "content": [{"type": "text", "text": content}]}
        if isinstance(content, dict) and "type" in content: # Single content block
            return {"role": role, "content": [content]}
        if isinstance(content, list): # Already a list of content blocks
            return {"role": role, "content": content}
        raise ValueError(f"Unsupported content format for format_message: {type(content)}")

    def run_conversation(self,
                        user_message: str,
                        logger: logging.Logger,
                        history: Optional[List[Dict[str, Any]]] = None,
                        history_with_tools: Optional[List[Dict[str, Any]]] = None,
                        max_tokens: Optional[int] = None,
                        temperature: Optional[float] = None,
                        role1: str = 'user',
                        role2: str = 'assistant',
                        max_iterations: int = 10,
                        run_id: str = '',
                        run_iteration: int = 0):
        """
        Runs a conversation with the LLM, keeping 'history' minimal (text-only) and 'history_with_tools' full (with tools).
        
        Args:
            user_message: The user's input as a string
            logger: Logging instance
            history: Simplified conversation history (text-only), updated and returned
            history_with_tools: Full history including tool calls/results, used for LLM
            max_tokens: Maximum tokens for response
            temperature: Temperature for response
            role1: Role for user messages (default: 'user')
            role2: Role for assistant messages (default: 'assistant')
            max_iterations: Maximum iterations for tool call loops
            run_id: Unique run identifier
        
        Returns:
            ClaudeResponse with simplified history, content, and tool results
        """
        # Initialize histories if not provided
        if history is None:
            history = []

        user_msg_formatted = self.format_message(role1, user_message)
        history.append(user_msg_formatted)

        if history_with_tools is None:
            history_with_tools = copy.deepcopy(history)

        # Set parameters
        temperature = temperature if temperature is not None else self.temperature
        max_tokens = max_tokens if max_tokens is not None else self.max_tokens

        total_price = 0.0
        results_list = []

        for _ in range(max_iterations):
            # Get LLM response using full history
            response, price = self._get_completion(history_with_tools, temperature, max_tokens, logger)
            logger.info(f'response {response}')
            total_price += price or 0.0
            print(f'response {response.choices}')
            
            # Check if response has no choices (e.g., content was blocked by safety filters)
            if not response.choices:
                error_msg = "No response choices returned - content may have been blocked by safety filters"
                logger.warning(error_msg)
                history.append(self.format_message(role2, error_msg))
                return ClaudeResponse(
                    content=error_msg,
                    content_blocks=[{"type": "text", "text": error_msg}],
                    history=history,
                    results=results_list,
                    price=total_price
                )
            
            llm_message = response.choices[0].message

            # Extract text content for simplified history and ClaudeResponse
            text_content = ""
            content_blocks = []
            if llm_message.content:
                if isinstance(llm_message.content, str):
                    text_content = llm_message.content
                    content_blocks.append({"type": "text", "text": text_content})
                elif isinstance(llm_message.content, list):
                    content_blocks.extend(llm_message.content)
                    for block in llm_message.content:
                        if block.get("type") == "text" and not text_content:
                            text_content = block["text"]

            # Add assistant message to history_with_tools (full details)
            assistant_msg = {"role": role2, "content": llm_message.content}
            if llm_message.tool_calls:
                assistant_msg["tool_calls"] = [tc.model_dump() for tc in llm_message.tool_calls]
            history_with_tools.append(assistant_msg)

            # Handle tool calls
            if llm_message.tool_calls:
                # Add text content to simplified history (if any, e.g., Claude-style)
                if text_content:
                    history.append(self.format_message(role2, text_content))

                # Process each tool call
                for tool_call in llm_message.tool_calls:
                    tool_name = tool_call.function.name
                    tool_id = tool_call.id
                    tool_args = json.loads(tool_call.function.arguments) if isinstance(tool_call.function.arguments, str) else tool_call.function.arguments
                    try:
                        tool_block = type('ToolBlock', (), {'name': tool_name, 'input': tool_args})()
                        result_dict = self._execute_tool_calls(tool_block, run_id, run_iteration)
                        logger.info(f'result_dict {result_dict}')
                        result = result_dict.get(tool_name, {"error": f"Result missing for {tool_name}"})

                        # Add tool result to history_with_tools
                        history_with_tools.append({
                            "role": "tool",
                            "tool_call_id": tool_id,
                            "content": json.dumps(result) if not isinstance(result, str) else result
                        })
                        results_list.append({"function_call": {"type": "tool_use", "id": tool_id, "name": tool_name, "input": tool_args}, "result": result})
                    except Exception as e:
                        error_msg = f"Tool '{tool_name}' error: {str(e)}"
                        history_with_tools.append({"role": "tool", "tool_call_id": tool_id, "content": json.dumps({"error": error_msg})})
                        results_list.append({"function_call": {"type": "tool_use", "id": tool_id, "name": tool_name, "input": tool_args}, "error": error_msg})
                continue  # Loop back to process tool results

            # No tool calls: final response
            if text_content:
                history.append(self.format_message(role2, text_content))
            self.tool_history.extend(history_with_tools)
            logger.info(f'history_with_tools {history_with_tools}')
            return ClaudeResponse(
                content=text_content,
                content_blocks=content_blocks,
                history=history,
                results=results_list,
                price=total_price
            )

        # Max iterations reached
        error_msg = "Maximum iterations reached"
        history.append(self.format_message(role2, error_msg))
        return ClaudeResponse(
            content=error_msg,
            content_blocks=[{"type": "text", "text": error_msg}],
            history=history,
            results=results_list,
            price=total_price
        )

    def _execute_tool_calls(self, tool_block_singular: Any, run_id: str, run_iteration: int) -> dict:
        # tool_block_singular is the TempToolBlock with .name and .input
        results = {}
        tool_name = tool_block_singular.name
        arguments = tool_block_singular.input
        arguments["run_iteration"] = run_iteration
        arguments['agent_name'] = self.title

        logger = logging.getLogger(__name__)
        logger.info(f"Executing tool '{tool_name}' with arguments: {arguments}")
        try:
            # This should return the direct result of the tool, not nested under tool_name again
            tool_output = tool_registry.execute(tool_name, arguments, run_id)
            results[tool_name] = tool_output
        except Exception as e:
            logger.error(f"Error executing tool {tool_name}: {str(e)}")
            import traceback
            traceback.print_exc()
            results[tool_name] = {"error": str(e), "details": "Exception during tool_registry.execute"}
        return results # Returns a dict like {"tool_name": actual_tool_result}
    
    def extract_thinking(self, text: str) -> str:
        thinking_match = re.search(r'<thinking>(.*?)</thinking>', text, re.DOTALL)
        return thinking_match.group(1).strip() if thinking_match else None

    
    def extract_outside_thinking(self, text: str) -> str:
        thinking_pattern = r'<thinking>.*?</thinking>'
        outside_thinking = re.sub(thinking_pattern, '', text, flags=re.DOTALL)
        return outside_thinking.strip() if outside_thinking else None
    
    def extract_tag_content(self, text: str, tag: str) -> str:
        """Extract content between XML-like tags.
        
        Args:
            text: The text containing tags
            tag: The tag name to search for (without < > or / )
        
        Returns:
            The content between the opening and closing tags
        """
        pattern = f"<{tag}>(.*?)</{tag}>"
        match = re.search(pattern, text, re.DOTALL)  # re.DOTALL allows matching across newlines
        return match.group(1).strip() if match else ""

    def calculate_price(self, input_tokens: int, output_tokens: int) -> float:
        """Calculate the price for Claude 3 Sonnet usage.
        
        Args:
            input_tokens: Number of input tokens
            output_tokens: Number of output tokens
            
        Returns:
            float: Total price in USD
        """
        INPUT_PRICE_PER_MILLION = 3.0  # $3 per million tokens for input
        OUTPUT_PRICE_PER_MILLION = 15.0  # $15 per million tokens for output
        
        input_cost = (input_tokens / 1_000_000) * INPUT_PRICE_PER_MILLION
        output_cost = (output_tokens / 1_000_000) * OUTPUT_PRICE_PER_MILLION
        
        return input_cost + output_cost

    def create_image_content_block(self, image_source: Union[str, bytes], mime_type: Optional[str] = None) -> Dict[str, Any]:
        """
        Create an image content block for Claude API from either a URL or base64-encoded image.
        
        Args:
            image_source: Either a URL to an image or base64-encoded image data
            mime_type: MIME type of the image (required for base64 images)
            
        Returns:
            A formatted image content block
        """
        # Check if the source is a URL
        if isinstance(image_source, str) and (image_source.startswith('http://') or 
                                             image_source.startswith('https://')):
            return {
                "type": "image",
                "source": {
                    "type": "url",
                    "url": image_source
                }
            }
        # Handle base64-encoded image
        elif isinstance(image_source, str) and image_source.startswith('data:'):
            # Extract MIME type and base64 data from the data URI
            parts = image_source.split(';base64,')
            if len(parts) != 2:
                raise ValueError("Invalid data URI format")
            
            detected_mime = parts[0].replace('data:', '')
            base64_data = parts[1]
            
            return {
                "type": "image",
                "source": {
                    "type": "base64",
                    "media_type": mime_type or detected_mime,
                    "data": base64_data
                }
            }
        # Handle raw base64 string (no data URI prefix)
        elif isinstance(image_source, str):
            if not mime_type:
                raise ValueError("MIME type is required for raw base64 image data")
            
            return {
                "type": "image",
                "source": {
                    "type": "base64",
                    "media_type": mime_type,
                    "data": image_source
                }
            }
        # Handle bytes
        elif isinstance(image_source, bytes):
            if not mime_type:
                raise ValueError("MIME type is required for image bytes")
            
            base64_data = base64.b64encode(image_source).decode('utf-8')
            return {
                "type": "image",
                "source": {
                    "type": "base64",
                    "media_type": mime_type,
                    "data": base64_data
                }
            }
        else:
            raise ValueError(f"Unsupported image source type: {type(image_source)}")

    def create_multimodal_message(self, text: str, images: List[Union[str, bytes]], mime_types: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """
        Create a multimodal message containing text and images.
        
        Args:
            text: Text content
            images: List of image sources (URLs or base64/bytes)
            mime_types: Optional list of MIME types for base64/bytes images
            
        Returns:
            List of content blocks ready for Claude API
        """
        # Start with text block
        content_blocks = [{"type": "text", "text": text}]
        
        # Add image blocks
        for i, image in enumerate(images):
            mime_type = None if mime_types is None else mime_types[i] if i < len(mime_types) else None
            content_blocks.append(self.create_image_content_block(image, mime_type))
        
        return content_blocks

    def send_message_with_images(self, 
                               text: str, 
                               images: List[Union[str, bytes]], 
                               history: Optional[List[Dict[str, Any]]] = None,
                               mime_types: Optional[List[str]] = None,
                               **kwargs) -> ClaudeResponse:
        """
        Convenience method to send a message with text and images in one step.
        
        Args:
            text: Text content
            images: List of image sources (URLs or base64/bytes)
            history: Optional conversation history
            mime_types: Optional list of MIME types for base64/bytes images
            **kwargs: Additional arguments to pass to run_conversation
            
        Returns:
            Claude response with content and history
        """
        # Create multimodal message
        message_content = self.create_multimodal_message(text, images, mime_types)
        
        # Run conversation with this message
        return self.run_conversation(
            user_message=message_content,
            history=history,
            **kwargs
        )

