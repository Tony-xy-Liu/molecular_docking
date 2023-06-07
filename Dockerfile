FROM condaforge/mambaforge as build-env
# https://mamba.readthedocs.io/en/latest/user_guide/mamba.html

COPY ./lib/conda_env.yml /opt/
RUN mamba env create --no-default-packages -f /opt/conda_env.yml \
    && mamba clean -afy

# Singularity uses tini, but raises warnings
# we set it up here correctly for singularity
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

## We do some umask munging to avoid having to use chmod later on,
## as it is painfully slow on large directores in Docker.
RUN old_umask=`umask` && \
    umask 0000 && \
    umask $old_umask

# move to clean execution environment
# jammy is 22.04 long term support
# FROM ubuntu:jammy
FROM eclipse-temurin:17-jre-jammy
COPY ./lib /opt
COPY ./test_data /test_data
COPY ./src /app
COPY ./bin /app/bin
COPY --from=build-env /opt/conda/envs/for_container /opt/conda/envs/for_container
COPY --from=build-env /tini /tini
ENV PATH=/app/bin:/opt/conda/envs/for_container/bin:/opt/adfr/bin:/opt/p2rank:$PATH

# # openjdk-17-jre \
# RUN DEBIAN_FRONTEND=noninteractive \
#     && apt-get update && apt-get install -y --no-install-recommends \
#     default-jre-headless \
#     && apt-get clean \
#     && rm -rf /var/lib/apt/lists/*

# singularity doesn't use the -s flag, and that causes warnings
ENTRYPOINT ["/tini", "-s", "--"]
