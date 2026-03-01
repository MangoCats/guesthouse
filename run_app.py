"""Launch the ADU Editor web application.

Usage:
    python run_app.py [--port PORT] [--host HOST]

Opens the interactive building editor at http://localhost:5000
"""
import argparse
import os
import sys
import webbrowser
import threading

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _DIR)


def main():
    parser = argparse.ArgumentParser(description="ADU Editor")
    parser.add_argument("--port", type=int, default=5000)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--no-browser", action="store_true")
    args = parser.parse_args()

    from app.server import create_app
    app = create_app()

    if not args.no_browser:
        url = f"http://{args.host}:{args.port}"
        threading.Timer(1.5, lambda: webbrowser.open(url)).start()
        print(f"\n  ADU Editor starting at {url}\n")

    app.run(host=args.host, port=args.port, debug=True, use_reloader=False, threaded=True)


if __name__ == "__main__":
    main()
