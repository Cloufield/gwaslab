import re
SNPID_PATTERN = r'^\w+[:_-]\d+[:_-][ATCG]+[:_-][ATCG]+$'

#SNPID_PATTERN_EXTRACT = r'^(?:chr)?(\w+)[:_-](\d+)[:_-]([ATCG]+)[:_-]([ATCG]+)$'
SNPID_PATTERN_EXTRACT = (
    r'^(?P<CHR_PREFIX>chr)?'
    r'(?P<CHR>\w+)[:_-]'
    r'(?P<POS>\d+)[:_-]'
    r'(?P<NEA>[ACGT]+)[:_-]'
    r'(?P<EA>[ACGT]+)$'
)

SNPID_PATTERN_STRIP =(
    r'^(?P<CHR_PREFIX>chr)?'
    r'(?P<CHR>\w+)'
    r'[:_-]'
    r'(?P<POS>\d+)'
    r'[:_-]'
    r'(?P<REF>[ACGT]+)'
    r'[:_-]'
    r'(?P<ALT>[ACGT]+)$'
)

RSID_PATTERN = r'^rs\d+$'

CHR_PATTERN_EXTRACT = r'^(chr)?(\d{1,3}|[XYM]|MT)$'

SNPID_PATTERN_CONTAINS = r'[:_-]?\w+[:_-]\d+[:_-][ATCG]+[:_-][ATCG]+[:_-]?'


FLAGS=re.IGNORECASE|re.ASCII
CHR_PREFIX_PATTERN = r'^chr'
SNPID_SEP_PATTERN = r'[_-]'
