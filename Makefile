coverage:
	nosetests tests/test_functions.py --with-coverage
	coverage report
	coverage html
