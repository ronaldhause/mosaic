# mosaic pipeline container.
# base image pins R 4.3.2 for reproducibility.
FROM rocker/r-ver:4.3.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        samtools \
        bedtools \
        perl \
        wget \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# R deps. keep this list minimal — only what the scripts actually import.
RUN R -e "install.packages(c('data.table','optparse'), repos='https://cloud.r-project.org')"

# python deps for the build_reference.nf and msings-comparator scripts.
RUN pip3 install --no-cache-dir pandas

# misa. the perl script is distributed freely from
# https://webblast.ipk-gatersleben.de/misa/. we mirror by the published
# tarball; swap the URL for a pinned release if the default stops working.
RUN mkdir -p /opt/misa \
 && wget -q -O /opt/misa/misa.pl \
      https://webblast.ipk-gatersleben.de/misa/download/misa.pl \
 && chmod +x /opt/misa/misa.pl
ENV PATH="/opt/misa:${PATH}"

# mSINGS. the pypi package name may differ across forks; if `pip install msings`
# fails, swap in a git install, e.g.:
#   pip install git+https://bitbucket.org/uwlabmed/msings.git
RUN pip3 install --no-cache-dir msings || \
    echo "WARNING: msings not available on pypi — install manually in a derived image"

# annovar is not installed here (personal registration required). mount an
# install at runtime:
#   -v /path/to/annovar:/opt/annovar
# and pass --annovar_dir /opt/annovar --annovar_db /opt/annovar/humandb to
# build_reference.nf.

COPY bin/ /opt/mosaic/bin/
RUN chmod +x /opt/mosaic/bin/*.R /opt/mosaic/bin/*.py 2>/dev/null || true
ENV PATH="/opt/mosaic/bin:${PATH}"

WORKDIR /work
