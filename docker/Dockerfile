FROM python:3.9
RUN apt update

# install tabix
RUN apt install -y tabix

# install fonts used for plotting
RUN apt install -y cabextract libfontenc1 libmspack0 xfonts-encodings xfonts-utils
RUN wget http://ftp.de.debian.org/debian/pool/contrib/m/msttcorefonts/ttf-mscorefonts-installer_3.8_all.deb
RUN dpkg -i ttf-mscorefonts-installer_3.8_all.deb
RUN rm ttf-mscorefonts-installer_3.8_all.deb

# install gwaslab and jupyter notebook
RUN pip install --upgrade pip
RUN python -m pip install "gwaslab==3.4.39"
RUN python -m pip install notebook
