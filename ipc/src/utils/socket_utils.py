from typing import Any, Callable, TypeVar
from socket import socket, AF_UNIX, SOCK_STREAM
from contextlib import contextmanager
import json

T = TypeVar('T')

def send_request(sock: socket, request: dict) -> Any:
    """
    Send request to server and receive response

    Args:
        sock: Connected socket
        request: Request data to send

    Returns:
        Parsed JSON response

    Raises:
        RuntimeError: If response is empty
    """
    request_json = json.dumps(request)
    request_bytes = request_json.encode()
    request_length = len(request_bytes)
    print(f"Sending request JSON: {request_json}")
    print(f"Request size: {request_length} bytes")
    
    try:
        # Send 4-byte length prefix in big-endian format
        length_bytes = request_length.to_bytes(4, byteorder='big')
        sock.sendall(length_bytes)
        
        # Send the actual request data
        sock.sendall(request_bytes)
        print("Request sent, waiting for response...")

        # Read the 4-byte length prefix
        length_bytes = sock.recv(4)
        if len(length_bytes) < 4:
            raise RuntimeError("Failed to read response length")
        length = int.from_bytes(length_bytes, 'big')
        print(f"Response length: {length} bytes")

        # Read exactly 'length' bytes
        response = b""
        while len(response) < length:
            chunk = sock.recv(min(4096, length - len(response)))
            if not chunk:
                raise RuntimeError("Connection closed before full response")
            response += chunk
            print(f"Received {len(response)}/{length} bytes")

        print(f"Complete response received ({len(response)} bytes)")
        return json.loads(response.decode())
    except Exception as e:
        print(f"Error during request/response: {e}")
        print(f"Last known socket state: {sock}")
        raise

@contextmanager
def connect_socket(socket_path: str = "/tmp/gurnemanz.sock") -> socket:
    """
    Context manager for socket connections

    Args:
        socket_path: Path to Unix domain socket

    Yields:
        Connected socket object
    """
    sock = socket(AF_UNIX, SOCK_STREAM)
    try:
        sock.connect(socket_path)
        yield sock
    finally:
        sock.close()

def with_connection(socket_path: str, operation: Callable[[socket], T]) -> T:
    """
    Execute operation with a socket connection

    Args:
        socket_path: Path to Unix domain socket
        operation: Function to execute with socket

    Returns:
        Result of operation
    """
    with connect_socket(socket_path) as sock:
        return operation(sock)
