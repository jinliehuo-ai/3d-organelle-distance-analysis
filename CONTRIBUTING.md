# Contributing

Do not open pull requests containing raw microscopy data, result tables, microscope metadata, private paths, or participant information.

For algorithm changes:

1. Add a versioned script or configuration.
2. Describe the changed assumption and its numerical consequences.
3. Add or update a small synthetic test.
4. Run `python scripts/qa/check_public_release.py .`.
5. Include the exact software versions used for validation.
