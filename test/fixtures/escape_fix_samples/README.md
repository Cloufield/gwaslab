# Escape-sequence fix samples

Minimal before/after examples copied from gwaslab patterns. Used to validate
fix approaches **before** editing production modules.

Run:

```bash
# Compile samples with SyntaxWarning as error
python -W error::SyntaxWarning -m compileall -q test/fixtures/escape_fix_samples

# Behavioral + compile checks
pytest test/test_escape_sequence_fix_patterns.py -v
```

Files:

| File | Purpose |
|------|---------|
| `bad_patterns.py.txt` | Original offending literals (reference only; `.txt` so compileall skips) |
| `good_patterns.py` | Fixed equivalents; must compile cleanly under `-W error::SyntaxWarning` |
| `demo_all_patterns.py` | Runnable script printing runtime verification |
