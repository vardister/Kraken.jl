# Formatting gate for CI (.github/workflows/Format.yml). Also the way to format by hand, from the
# repo root:
#
#     julia .github/format_check.jl
#
# Note it formats **in place** and then reports — that is how JuliaFormatter works, and it is why
# this doubles as the fix command. Settings come from .JuliaFormatter.toml (blue, indent 4,
# margin 120).
#
# Explicit paths rather than ".": JuliaFormatter does not read .gitignore, so `format(".")` descends
# into .git/ and tries to parse ref files — this repo has branches whose names end in `.jl`, so that
# produces a screenful of bogus ParseError warnings.
using JuliaFormatter

const PATHS = ["src", "ext", "test", "dev", "docs", "examples", "benchmark"]

targets = filter(isdir, PATHS)
append!(targets, filter(f -> endswith(f, ".jl"), readdir(".")))

if all(format, targets)
    println("All files are formatted.")
else
    # Show what would change, so the failure says what to fix rather than only that something is.
    write(stdout, read(`git diff`, String))
    @error "Some files are not formatted. Re-run `julia .github/format_check.jl` and commit the result."
    exit(1)
end
