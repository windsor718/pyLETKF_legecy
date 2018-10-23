#!/opt/local/bin/python
# -*- coding: utf-8 -*-
import datetime
import time
import requests
import json
import sys

class sendSlack(object):

    def __init__(self,url):

        self.webhookUrl = url


    def success(self,user,text):

        now   = datetime.datetime.now()
        posix = int(time.mktime(now.timetuple()))
        requests.post(self.webhookUrl, data = json.dumps({
            "attachments": [
            {
                "fallback": "New notification from %s [green]" % user,
                "color": "#27AE60",
                "fields": [
                    {
                        "title": "Success",
                        "value": text,
                        "short": "false",
                    }
                ],
                "ts": posix,
            }
        ],
        "username": user,
        "link_name": 1
        }))


    def failed(user,text):
        now   = datetime.datetime.now()
        posix = int(time.mktime(now.timetuple()))
        requests.post(self.webhookUrl, data = json.dumps({
            "attachments": [
            {
                "fallback": "New notification from %s [red]" % user,
                "color": "#C70039",
                "fields": [
                    {
                        "title": "Error",
                        "value": text,
                        "short": "false",
                    }
                ],
                "ts": posix,
            }
        ],
        "username": user,
        "link_name": 1
        }))


    def progress(user,text):
        now   = datetime.datetime.now()
        posix = int(time.mktime(now.timetuple()))
        requests.post(self.webhookUrl, data = json.dumps({
            "attachments": [
            {
                "fallback": "New notification from %s [blue]" % user,
                "color": "#3498DB",
                "fields": [
                    {
                        "title": "InProgress",
                        "value": text,
                        "short": "false",
                    }
                ],
                "ts": posix,
            }
        ],
        "username": user,
        "link_name": 1
        }))


    def unknown(user,text):
        now   = datetime.datetime.now()
        posix = int(time.mktime(now.timetuple()))
        requests.post(self.webhookUrl, data = json.dumps({
            "attachments": [
            {
                "fallback": "New notification from %s [unknown]" % user,
                "pretext": "I got unknown flag[%s]. Displayed as unknown." % flag,
                "color": "#34495E",
                "fields": [
                    {
                        "title": "Unknown",
                        "value": text,
                        "short": "false",
                    }
                ],
                "ts": posix,
            }
        ],
        "username": user,
        "link_name": 1
        }))
