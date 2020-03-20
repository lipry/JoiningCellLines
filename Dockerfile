FROM python:3

RUN addgroup --gid 1000 stud
RUN adduser --disabled-password --gecos "" --force-badname --gid 1000 lipreri
USER lipreri

ADD ./ app/
COPY requirements.txt app/
WORKDIR app

RUN pip install --no-cache-dir -r requirements.txt

CMD [ "python", "./main.py" ]
