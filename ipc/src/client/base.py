from typing import Any, Callable, TypeVar
from socket import socket, AF_UNIX, SOCK_STREAM
from contextlib import contextmanager
from ..utils.socket_utils import send_request

T = TypeVar('T')

@contextmanager
def connect_socket(socket_path: str = "/tmp/gurnemanz.sock") -> socket:
    """Context manager for socket connection"""
    sock = socket(AF_UNIX, SOCK_STREAM)
    try:
        sock.connect(socket_path)
        yield sock
    finally:
        sock.close()

def with_connection(socket_path: str, operation: Callable[[socket], T]) -> T:
    """Helper to run an operation with a socket connection"""
    with connect_socket(socket_path) as sock:
        return operation(sock)
