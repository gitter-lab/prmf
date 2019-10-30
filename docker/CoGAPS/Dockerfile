
FROM rocker/r-base:latest

# Install CoGAPS and its dependencies
RUN useradd -m CoGAPS_user
COPY ./CoGAPS_setup.R /home/CoGAPS_user/
COPY ./CoGAPS_wrapper.R /home/CoGAPS_user/
WORKDIR /home/CoGAPS_user/
RUN ["Rscript", "CoGAPS_setup.R"]
ENV PATH="/home/CoGAPS_user:${PATH}"
ENTRYPOINT ["/home/CoGAPS_user/CoGAPS_wrapper.R", "--help"]
