

import unittest

from traits.api import Any, HasTraits, Unicode

from ..splash_screen_log_handler import SplashScreenLogHandler


class DummySplashScreen(HasTraits):
    text = Unicode('original')


class DummyRecord(object):
    def __init__(self, message):
        self.message = message

    def getMessage(self):
        return self.message


class TestSplashScreenLogHandler(unittest.TestCase):

    def setUp(self):
        self.ss = DummySplashScreen()
        self.sslh = SplashScreenLogHandler(self.ss)

    def test_unicode_message(self):
        self.assertEqual(self.ss.text, 'original')
        message = 'G\u00f6khan'
        self.sslh.emit(DummyRecord(message))
        self.assertEqual(self.ss.text, message + '...')

    def test_ascii_message(self):
        message = 'Goekhan'
        self.sslh.emit(DummyRecord(message))
        self.assertEqual(self.ss.text, message + '...')
