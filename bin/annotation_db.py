"""Connection-file validation shared by LEAP annotation tasks and preflight."""

import argparse
import json
import re
from pathlib import Path
import sys
from urllib.parse import quote, quote_plus


def load_connection(path):
    try:
        with Path(path).open(encoding='utf-8') as handle:
            value = json.load(handle)
    except (OSError, ValueError, UnicodeError):
        raise ValueError('Cannot read a valid database secrets JSON file') from None
    if not isinstance(value, dict):
        raise ValueError('Database secrets JSON must contain an object')
    for name in ('host', 'user'):
        if not isinstance(value.get(name), str) or not value[name].strip():
            raise ValueError('Database secrets JSON requires a nonempty ' + name)
    database = value.get('database') or value.get('dbname')
    if not isinstance(database, str) or not database.strip():
        raise ValueError('Database secrets JSON requires database or dbname')
    port = value.get('port')
    if isinstance(port, bool) or not isinstance(port, (int, str)) or not str(port).isdigit() or not 1 <= int(port) <= 65535:
        raise ValueError('Database secrets JSON requires port in 1..65535')
    password = value.get('password', value.get('pass', ''))
    if not isinstance(password, str):
        raise ValueError('Database secrets password/pass must be a string (empty is allowed)')
    if value.get('dialect', 'mysql') != 'mysql':
        raise ValueError('LEAP annotation export requires the mysql dialect')
    for first, second in (('database', 'dbname'), ('password', 'pass')):
        if first in value and second in value and value[first] != value[second]:
            raise ValueError('Database secrets JSON has conflicting ' + first + '/' + second)
    return dict(host=value['host'], port=int(port), user=value['user'],
                password=password, database=database, dialect='mysql')


def redact_connection_text(message, connection):
    """Retain useful schema/sequence errors without exposing connection values."""
    values = set()
    for value in connection.values():
        if value is not None and str(value):
            values.update((str(value), quote(str(value), safe=''), quote_plus(str(value))))
    for value in sorted(values, key=len, reverse=True):
        message = re.sub(re.escape(value), '<redacted>', message, flags=re.IGNORECASE)
    return message[:6000]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('db_secrets')
    args = parser.parse_args()
    try:
        load_connection(args.db_secrets)
        return 0
    except ValueError as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
