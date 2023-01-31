install:
	pip install --upgrade pip
	pip install -r requirements.txt
format:
	black *.py fragment_analyzer/*.py
lint:
	pylint --disable=R,C *.py fragment_analyzer/*.py