NAME=molecular_docking
DOCKER_IMAGE=quay.io/hallamlab/$NAME
echo image: $DOCKER_IMAGE
echo ""

HERE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

case $1 in
    --build|-b)
        # change the url in python if not txyliu
        # build the docker container locally *with the cog db* (see above)
        docker build -t $DOCKER_IMAGE .
    ;;
    --push|-p)
        # login and push image to quay.io, remember to change the python constants in src/
        # sudo docker login quay.io
	    docker push $DOCKER_IMAGE:latest
    ;;
    --sif)
        # test build singularity
        singularity build $NAME.sif docker-daemon://$DOCKER_IMAGE:latest
    ;;
    # --pip-setup)
    #     # make an environment before hand
    #     # in that env, install these build tools
    #     pip install build
    # ;;
    # --pip-build|-l)
    #     # build the packge for container
    #     # this is NOT suitable for pypi because of poorly handled dependencies
    #     rm -r build && rm -r dist
    #     python -m build --wheel
    # ;;
    # --pip-remove|-x)
    #     pip uninstall $NAME -y
    # ;;
    --run|-r)
        # test run docker image
            # --mount type=bind,source="$HERE/scratch",target="/ws" \
            # --mount type=bind,source="$HERE/scratch/res",target="/ref"\
            # --mount type=bind,source="$HERE/scratch/res/.ncbi",target="/.ncbi" \
            # --mount type=bind,source="$HERE/test",target="/ws" \
            # --mount type=bind,source="$HERE/test/checkm_db",target="/checkm_db" \
            # -e XDG_CACHE_HOME="/ws"\
            
        docker run -it --rm \
            --mount type=bind,source="$HERE",target="/ws" \
            --workdir="/ws" \
            -u $(id -u):$(id -g) \
            $DOCKER_IMAGE \
            /bin/bash 
    ;;
    -t)
        out=cache/MD_test_results
        rm -r $out
            # -e XDG_CACHE_HOME="/ws"\
            # moldo --attempts 10 --active-sites 1 --flexes 0 --exhaustiveness 8 --threads 14 --seed 12345 \
        docker run -it --rm \
            --mount type=bind,source="$HERE/src/moldo",target="/opt/conda/envs/for_container/lib/python3.9/site-packages/moldo" \
            --mount type=bind,source="$HERE/lib/adfr",target="/opt/adfr" \
            --mount type=bind,source="$HERE/lib/p2rank",target="/opt/p2rank" \
            --mount type=bind,source="$HERE",target="/ws" \
            --workdir="/ws" \
            -u $(id -u):$(id -g) \
            $DOCKER_IMAGE \
            moldo --attempts 10 --active-sites 1 --flexes 4 --exhaustiveness 16 --threads 14 --seed 12345 \
            -r /test_data/malate_dehydrogenase/mdh_3hhp.pdb -l /test_data/malate_dehydrogenase/s_malate_ChEBI_15589.sdf \
            -o $out
    ;;
    *)
        echo "bad option"
        echo $1
    ;;
esac
