#!/bin/bash
python3 manage.py runserver 7000 &
python3 manage.py process_tasks
