JL := "julia --project"

# Default recipe
default: init test

# Initialize the project
init:
    {{JL}} -e 'using Pkg; Pkg.precompile()'

# Initialize docs environment
init-docs:
    {{JL}} -e 'using Pkg; Pkg.activate("docs"); Pkg.precompile()'

# Update packages
update:
    {{JL}} -e 'using Pkg; Pkg.update(); Pkg.precompile()'

# Update packages in docs environment
update-docs:
    {{JL}} -e 'using Pkg; Pkg.activate("docs"); Pkg.update(); Pkg.precompile()'

# Run tests
test:
    {{JL}} -e 'using Pkg; Pkg.test()'

# Serve documentation
servedocs:
    {{JL}} -e 'using Pkg; Pkg.activate("docs"); using LiveServer; servedocs(;skip_dirs = ["docs/src/assets", "docs/src/generated"])'

# Clean build artifacts
clean:
    rm -rf docs/build
    find . -name "*.cov" -type f -print0 | xargs -0 /bin/rm -f
