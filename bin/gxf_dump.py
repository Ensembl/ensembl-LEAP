#!/usr/bin/env python3
"""Bundled db2gxf CLI with credential-safe live-export diagnostics."""

from contextlib import redirect_stderr, redirect_stdout
import io
import sys

sys.dont_write_bytecode = True

from gxf_dump.exporter import export_main
from annotation_db import load_connection, redact_connection_text


def main():
    if not any(arg.split('=', 1)[0] == '--db-secrets' for arg in sys.argv[1:]):
        return export_main()
    output, errors = io.StringIO(), io.StringIO()
    # Validate separately so malformed JSON cannot leak a source fragment through
    # a decoder traceback. Subsequent driver diagnostics can be redacted safely.
    try:
        arguments = sys.argv[1:]
        option = next(i for i, value in enumerate(arguments) if value.split('=', 1)[0] == '--db-secrets')
        filename = arguments[option].split('=', 1)[1] if '=' in arguments[option] else arguments[option + 1]
        connection = load_connection(filename)
    except (ValueError, IndexError):
        print('ERROR: Invalid database secrets JSON; check the connection file and required fields.', file=sys.stderr)
        return 2
    try:
        with redirect_stdout(output), redirect_stderr(errors):
            result = export_main()
    except (Exception, SystemExit) as error:
        if isinstance(error, SystemExit) and error.code == 0:
            print(output.getvalue(), end='')
            return 0
        details = errors.getvalue().strip() if isinstance(error, SystemExit) else str(error)
        print('ERROR: Database export failed: ' + redact_connection_text(details, connection), file=sys.stderr)
        return 2
    print(output.getvalue(), end='')
    if errors.getvalue().strip():
        print(redact_connection_text(errors.getvalue().strip(), connection), file=sys.stderr)
    return result


if __name__ == '__main__':
    sys.exit(main())
