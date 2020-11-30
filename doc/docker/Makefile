all : thinkbox

thinkbox:
	docker build -t dxcompiler_workbench .
	docker run \
           --name dxcompiler_thinkbox \
           --privileged \
	   --volume /tmp:/tmp \
	   --volume /var/run/docker.sock:/var/run/docker.sock \
	-it dxcompiler_workbench

stop:
	docker stop dxcompiler_thinkbox || true
	docker rm dxcompiler_thinkbox || true
