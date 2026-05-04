.PHONY: help install run-scenarios run-validation run-spatial-model run-monte-carlo run-results run-tests docker-up docker-down clean

help:
	@echo "Chernobyl RBMK-1000 Accident Simulation — Build & Run Commands"
	@echo "=============================================================="
	@echo ""
	@echo "Scenarios & Analysis:"
	@echo "  make run-scenarios         Run all 9 scenarios in Julia"
	@echo "  make run-monte-carlo       Run 500-run Monte Carlo UQ"
	@echo "  make run-validation        Run SKALA data validation"
	@echo "  make run-spatial-model     Run 2D spatial PDE model (long, ~30min)"
	@echo "  make run-results           Generate full results package"
	@echo "  make run-tests             Run scenario smoke tests"
	@echo ""
	@echo "Docker Deployment:"
	@echo "  make docker-up             Start all services (docker-compose)"
	@echo "  make docker-down           Stop all services"
	@echo ""
	@echo "Utilities:"
	@echo "  make install               Install all dependencies"
	@echo "  make clean                 Clean generated files"
	@echo ""

install:
	@echo "Installing Julia dependencies (Project.toml)..."
	julia --project=. -e "using Pkg; Pkg.instantiate()"

run-scenarios:
	julia --project=. scripts/run_scenarios.jl

run-monte-carlo:
	julia --project=. scripts/run_monte_carlo.jl

run-validation:
	julia --project=. scripts/run_validation.jl

run-spatial-model:
	@echo "WARNING: 2D spatial model runs ~30-60 minutes for full solution"
	julia --project=. scripts/run_spatial_model.jl

run-results:
	julia --project=. scripts/run_results.jl

run-tests:
	julia --project=. tests/test_scenarios_smoke.jl
	julia --project=. tests/test_scenarios_compare.jl
	julia --project=. tests/test_scenarios_f.jl

docker-up:
	docker-compose up -d
	@echo "Services starting..."
	@sleep 5
	docker-compose logs -f

docker-down:
	docker-compose down

clean:
	rm -f *.csv *.png *.gif
	rm -rf __pycache__ .pytest_cache
	find . -name "*.jl~" -delete
	@echo "Cleaned up generated files"

.DEFAULT_GOAL := help
