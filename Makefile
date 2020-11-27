test: 1c.mdv
	python3 cover_by_cycles.py

1c.mdv:
	genice 1c -r 4 4 4 -f mdview > $@

prepare:
	pip install genice
