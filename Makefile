install:
	pip install --upgrade pip
	pip install -r requirements.txt
format:
	black *.py fragment_analyzer/*.py fragment_analyzer/ladders/* fragment_analyzer/reports/*.py
lint:
	pylint --disable=R,C *.py fragment_analyzer/*.py fragment_analyzer/reports/*.py
clean:
	rm -rf dist/ build/ *.egg-info
build:
	python -m build
