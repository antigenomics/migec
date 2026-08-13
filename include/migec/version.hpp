#ifndef MIGEC_VERSION_HPP
#define MIGEC_VERSION_HPP

// Kept in sync by hand with pyproject.toml and python/migec/__init__.py. The release checklist in
// CLAUDE.md lists all three; CI prints them so a mismatch is visible before a release, and
// publish.yml asserts the pyproject version equals the release tag.
#define MIGEC_VERSION "2.0.0a4"

#endif  // MIGEC_VERSION_HPP
