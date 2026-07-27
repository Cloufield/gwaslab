"""Shared gzip compression defaults for io writers (GNU ``gzip`` default level)."""

GZIP_COMPRESSLEVEL = 6

PANDAS_GZIP_COMPRESSION = {"method": "gzip", "compresslevel": GZIP_COMPRESSLEVEL}
