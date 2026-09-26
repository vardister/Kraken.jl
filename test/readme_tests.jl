using Test
using Markdown

# Every ```julia block in README.md, run in order in one fresh module, the way a reader would paste
# them into a REPL. Plan task 8.3: the old README's group-speed example called `find_kr` with a
# `method=Roots.A42()` keyword long after the solver had moved to NonlinearSolve, and nothing noticed.
#
# Only the language tag decides what runs. The install line is fenced `julia-repl` and the shell
# commands `sh`, so neither is executed; anything that needs a package the test environment lacks
# (Makie, say) belongs in prose, not in a `julia` block. The blocks are not doctests with checked
# output: GitHub renders a `jldoctest` fence as plain text, and a README is for reading on GitHub. The
# check is that they still run, which is the kind of drift that happened.

function readme_julia_blocks(path)
    blocks = String[]
    function collect!(node)
        if node isa Markdown.Code
            node.language == "julia" && push!(blocks, node.code)
        elseif hasproperty(node, :content)
            foreach(collect!, node.content)
        end
    end
    collect!(Markdown.parse_file(path))
    return blocks
end

@testset "README code blocks run" begin
    blocks = readme_julia_blocks(joinpath(@__DIR__, "..", "README.md"))
    # A parser change that silently found nothing would make this testset vacuous.
    @test length(blocks) >= 5

    sandbox = Module(:ReadmeSandbox)
    for (i, code) in enumerate(blocks)
        first_line = strip(first(split(code, '\n')))
        @testset "block $i: $first_line" begin
            @test (Base.include_string(sandbox, code, "README.md block $i"); true)
        end
    end
end
