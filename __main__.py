import sys

from .cli import main

if __name__ == "__main__":
    try:
        code = main()
    except BrokenPipeError:  # e.g. output piped into `head`
        code = 0
    sys.exit(code)
