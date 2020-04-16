
# Try to import the local version (because if it exists, it will be more recent than miri.apt_parser
try:
    import apt_parser
except ModuleNotFoundError:
    from miri import apt_parser

apt_parser.init_log()

observations = apt_parser.parse_apt("1028.aptx")

print(observations)
