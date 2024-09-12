# Definitions
repository = gcr.io/ucsd-medicine-cast
identifier = trtools-annotatr
version = 0.0.0
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = ${repository}/${identifier}
tag = ${version}--${git_commit}

# Steps
build:
	# do the docker build
	docker build --platform="linux/amd64" -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest

push:
	# Requires ~/.dockercfg
	docker push ${name}:${tag}