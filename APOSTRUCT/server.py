"""Standalone localhost server for the APOSTRUCT frontend and API."""

from __future__ import annotations

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import threading
from typing import Any
from urllib.parse import parse_qs, urlparse
import webbrowser

from APOSTRUCT import (
    _saved_debug_cases_available,
    handle_api,
    handle_api_post,
)


FRONTEND_ROOT = Path(__file__).resolve().parent / "Frontend"
CONTENT_TYPES = {
    ".csv": "text/csv; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".svg": "image/svg+xml",
}
MAX_REQUEST_BYTES = 16 * 1024 * 1024


def _asset_path(path: str) -> Path | None:
    if path in {"/", "/index.html", "/tools/APOSTRUCT", "/tools/APOSTRUCT/"}:
        relative = "index.html"
    elif path.startswith("/tools/APOSTRUCT/"):
        relative = path.removeprefix("/tools/APOSTRUCT/")
    else:
        relative = path.removeprefix("/")
    try:
        resolved = (FRONTEND_ROOT / relative).resolve()
        resolved.relative_to(FRONTEND_ROOT.resolve())
    except (OSError, ValueError):
        return None
    return resolved if resolved.is_file() else None


def make_handler():
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *_args: Any) -> None:
            return

        def _send(
            self,
            status: int,
            body: bytes | str | dict[str, Any] | list[Any],
            content_type: str = "text/plain; charset=utf-8",
        ) -> None:
            if isinstance(body, (dict, list)):
                body = json.dumps(body, ensure_ascii=False, allow_nan=False)
                content_type = "application/json; charset=utf-8"
            data = body.encode() if isinstance(body, str) else body
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(data)

        def _error(self, status: int, exc: Exception) -> None:
            self._send(status, {"error": f"{type(exc).__name__}: {exc}"})

        def do_GET(self) -> None:
            parsed = urlparse(self.path)
            if parsed.path.startswith("/api/APOSTRUCT/"):
                endpoint = parsed.path.removeprefix("/api/APOSTRUCT/")
                try:
                    self._send(200, handle_api(endpoint, parse_qs(parsed.query)))
                except (IndexError, KeyError, ValueError) as exc:
                    self._error(400, exc)
                return
            asset = _asset_path(parsed.path)
            if asset is None:
                self._send(404, "Not Found")
                return
            self._send(
                200,
                asset.read_bytes(),
                CONTENT_TYPES.get(asset.suffix.casefold(), "application/octet-stream"),
            )

        def do_POST(self) -> None:
            parsed = urlparse(self.path)
            if parsed.path == "/shutdown":
                self._send(200, "ok")
                threading.Thread(target=self.server.shutdown, daemon=True).start()
                return
            if not parsed.path.startswith("/api/APOSTRUCT/"):
                self._send(404, "Not Found")
                return
            try:
                length = int(self.headers.get("Content-Length", "0") or "0")
            except ValueError as exc:
                self._error(400, exc)
                return
            if not 0 <= length <= MAX_REQUEST_BYTES:
                self._send(413, {"error": "request body is too large"})
                return
            try:
                raw = self.rfile.read(length) if length else b"{}"
                payload = json.loads(raw.decode())
                if not isinstance(payload, dict):
                    raise ValueError("request body must be a JSON object")
                endpoint = parsed.path.removeprefix("/api/APOSTRUCT/")
                result = handle_api_post(endpoint, parse_qs(parsed.query), payload)
                self._send(200, result)
            except (
                AttributeError,
                IndexError,
                KeyError,
                TypeError,
                UnicodeDecodeError,
                ValueError,
                json.JSONDecodeError,
            ) as exc:
                self._error(400, exc)

    return Handler


def create_server(host: str, port: int) -> ThreadingHTTPServer:
    if not isinstance(port, int) or isinstance(port, bool) or not 0 <= port <= 65535:
        raise ValueError("port must be an integer from 0 through 65535")
    class Server(ThreadingHTTPServer):
        def server_close(self) -> None:
            if _saved_debug_cases_available():
                from APOSTRUCT.debug.jobs import cancel_all_debug_jobs

                cancel_all_debug_jobs()
            super().server_close()

    return Server((host, port), make_handler())


def serve(*, host: str = "127.0.0.1", port: int = 8300, open_browser: bool = False) -> None:
    server = create_server(host, port)
    actual_host, actual_port = server.server_address[:2]
    url = f"http://{actual_host}:{actual_port}/"
    print(f"APOSTRUCT: {url}")
    print("Stop: Ctrl+C or POST /shutdown")
    if open_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


__all__ = ["FRONTEND_ROOT", "create_server", "make_handler", "serve"]
