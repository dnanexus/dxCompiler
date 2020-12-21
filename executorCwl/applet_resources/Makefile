SHELL=/bin/bash -e

ifneq (,$(wildcard ./.env))
	include .env
	export
else
	# we can remove this once asset-builder supports .env files
	export NODEJS_VERSION=14
endif

all:
	# download nodejs installer
	curl -sL https://deb.nodesource.com/setup_${NODEJS_VERSION}.x | sudo -E bash -
	# install nodejs
	sudo apt-get install -y nodejs
	# remove the installer
	rm setup_14.x