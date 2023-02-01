install:
	pip install --upgrade pip
	pip install -r requirements.txt
format:
	black *.py fragment_analyzer/*.py fragment_analyzer/ladders/*
lint:
	pylint --disable=R,C *.py fragment_analyzer/*.py
