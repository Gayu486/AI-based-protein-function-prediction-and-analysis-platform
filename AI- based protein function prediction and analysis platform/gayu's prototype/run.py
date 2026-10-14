"""
Startup script for AI-Based Protein Function Prediction and Analysis Platform.
Runs the Flask application and opens the web browser automatically.
"""

import sys
import webbrowser
import threading
import time
from app import app

def open_browser(port):
    time.sleep(1.2)
    url = f"http://127.0.0.1:{port}"
    print(f"Opening browser at: {url}")
    webbrowser.open(url)

if __name__ == "__main__":
    port = 5000
    threading.Thread(target=open_browser, args=(port,), daemon=True).start()
    app.run(host="127.0.0.1", port=port, debug=False)
