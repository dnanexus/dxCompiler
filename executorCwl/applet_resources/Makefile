SHELL=/bin/bash -e

ifneq (,$(wildcard ./.env))
	include .env
	export
else
	# we can remove this once asset-builder supports .env files
	export NODEJS_VERSION=14.15.3
endif

PKG=node-v${NODEJS_VERSION}-linux-x64

all:
	# download nodejs
	echo "downloading ${PKG}"
	wget https://nodejs.org/dist/v${NODEJS_VERSION}/${PKG}.tar.xz
	# unpack nodejs tar
	tar -xf ${PKG}.tar.xz
	# move the executable
	mv ${PKG}/bin/node /usr/bin
	# remove the files
	rm -Rf ${PKG} ${PKG}.tar.xz