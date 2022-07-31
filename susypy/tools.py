import sys
import time
import functools

from typing import Any

def pretty_print(d):
	for key, value in d.items():
		if isinstance(value, dict):
			print(f'======= {key} =======')
			pretty_print(value)
		else:
			print(f'{key}: {value}')

def send_msg(account_sid: str, auth_token: str, from_: str, to: str, send_val: bool = False) -> None:
	"""Sends message when function ends"""

	def decorator_send_msg(func):
		@functools.wraps(func)

		def wrapper_send_msg(*args, **kwargs):
			from twilio.rest import Client

			client = Client(account_sid, auth_token)
			value = func(*args, **kwargs)
			meg_end = f'It returned {value}.' if send_val else 'GOOD LUCK :)'

			message = client.messages.create(
				body=f"SusyPy's {func.__name__!r} is done calculating! {meg_end}",
				from_=from_,
				to=to
			)

			print(message.sid)

			return value

		return wrapper_send_msg

	return decorator_send_msg

def timer(func) -> Any:
	"""Print the runtime of the decorated function"""

	@functools.wraps(func)

	def wrapper_timer(*args, **kwargs):
		start_time = time.perf_counter()
		value = func(*args, **kwargs)
		end_time = time.perf_counter()
		run_time = end_time - start_time

		print(f'Finished {func.__name__!r} in {run_time:.4f} secs')

		return value

	return wrapper_timer

def logger(file_name: str) -> Any:
	"""Saves print statments into a log file"""

	def decorator_logger(func):
		@functools.wraps(func)

		def wrapper_logger(*args, **kwargs):
			file_str = f'{file_name}_{round(time.time())}.txt'
			print(f'Saving to {file_str!r}...')
			temp = sys.stdout
			sys.stdout = open(file_str, 'a')
			value = func(*args, **kwargs)
			sys.stdout.close()
			sys.stdout = temp
			return value

		return wrapper_logger

	return decorator_logger