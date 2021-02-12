name = htp_md
tag = ${name}
build_tag= ${name}-build
container = ${build_tag}
workspace = $(shell pwd)
env = local

clean:
	rm -rf .coverage  coverage_reports

lint:
	pycodestyle htpmd
	flake8 --count --show-source --statistics htpmd
	flake8 --count --exit-zero --max-complexity=20 --statistics htpmd

test:
	pytest htpmd --color=yes --cov=htpmd --cov-config=.coveragerc --cov-report html:coverage_reports

clean-docker: clean
	docker rm -f ${name}  || true
	docker rm -f ${build_tag}  || true

build-docker: clean-docker
	docker build -t ${build_tag} .

test-docker: build-docker
	docker run -d --name ${container} ${build_tag}

	docker exec ${container} bash -c 'pycodestyle htpmd'
	docker exec ${container} bash -c 'flake8 --count --show-source --statistics htpmd'
	 # exit-zero treats all errors as warnings.
    docker exec ${container} bash -c 'flake8 --count --exit-zero --max-complexity=20 --statistics htpmd'
	docker exec ${container} bash -c 'pytest htpmd --color=yes --cov=htpmd --cov-config=.coveragerc --cov-report html:coverage_reports'

view-coverage:
	open ./coverage_reports/index.html 
