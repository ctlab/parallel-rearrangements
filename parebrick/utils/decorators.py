import time

sep = '-' * 80

def decorate(module_name, logger):
    def actual_decorator(func):

        def wrapper(*args, **kwargs):
            start = time.time()

            logger.info(f'<Running module: {module_name}>')
            print(sep)
            return_value = func(*args, **kwargs)

            end = time.time()

            print(sep)
            logger.info(f'<Module finished: {module_name}>, elapsed {end - start} seconds')
            print()
            print()
            return return_value

        return wrapper

    return actual_decorator

# def print_module_name():
