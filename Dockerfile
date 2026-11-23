# Ambuild's Poreblazer fork as a small runtime image: the executables at
# /opt/poreblazer/poreblazer.exe (serial, OpenMP) and /opt/poreblazer/poreblazer-mpi.exe
# (the grid split over MPI ranks; Open MPI), and the libraries they need, on Debian
# bookworm. Published to ghcr.io/st7ma784/poreblazer by .github/workflows/ci.yml once the
# tests pass. Other images take the executables from it:
#
#   COPY --from=ghcr.io/st7ma784/poreblazer:sha-<commit> /opt/poreblazer/ /opt/poreblazer/
#
# (poreblazer.exe needs libgfortran5 and libgomp1; poreblazer-mpi.exe also libopenmpi3,
# and openmpi-bin for mpirun). Or run it on a directory holding input.dat, defaults.dat,
# UFF.atoms and the structure:
#
#   docker run --rm -i -v "$PWD:/work" ghcr.io/st7ma784/poreblazer:ambuild < input.dat
#
# On a cluster, build poreblazer-mpi.exe against the site's MPI (make -C src mpi with its
# mpif90) rather than copying this one, so srun and the interconnect work.

FROM debian:bookworm-slim AS build
RUN apt-get update \
 && apt-get install -y --no-install-recommends gfortran make libopenmpi-dev \
 && rm -rf /var/lib/apt/lists/*
COPY src /src
# The Makefile's flags: -O2 -fopenmp -ffp-contract=off -fvect-cost-model=dynamic
RUN make -C /src && make -C /src mpi && strip /src/poreblazer.exe /src/poreblazer-mpi.exe

FROM debian:bookworm-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends libgfortran5 libgomp1 libopenmpi3 openmpi-bin \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/poreblazer.exe /src/poreblazer-mpi.exe /opt/poreblazer/
COPY FORK.md LICENSE /opt/poreblazer/
LABEL org.opencontainers.image.source="https://github.com/st7ma784/poreblazer" \
      org.opencontainers.image.description="Ambuild's Poreblazer fork: OpenMP and MPI, cell lists, exact percolation labelling; output identical to Poreblazer 3.0.5 by default" \
      org.opencontainers.image.licenses="GPL-3.0"
ENV PATH=/opt/poreblazer:$PATH
WORKDIR /work
ENTRYPOINT ["poreblazer.exe"]
