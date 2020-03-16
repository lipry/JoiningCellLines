FROM python:3

ADD ./ app/
COPY requirements.txt app/
WORKDIR app

RUN pip install --no-cache-dir -r requirements.txt

CMD [ "python", "./main.py" ]
