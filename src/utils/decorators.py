import time

sep = '-' * 80

def decorate(module_name):
    def actual_decorator(func):

        def wrapper(*args, **kwargs):
            start = time.time()

            print(f'<Running module: {module_name}>')
            print(sep, '')
            return_value = func(*args, **kwargs)

            end = time.time()

            print(sep, '')
            print(f'<Module finished: {module_name}>, elapsed {end - start} seconds')
            print()
            print()
            return return_value

        return wrapper

    return actual_decorator

# def print_module_name():
