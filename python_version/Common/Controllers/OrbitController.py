import json
import Utility.helpers as h
import numpy as np

class OrbitController:
    # Data instanced from Input Manifold
    _data = []
    _generation_id = None
    _table = 'orbit_generation'
    _parent_table = None
    _child_table = None

    # TODO Add SingletonMeta Here
    # Possibly add process ID path

    # Cache instanced from class loads
    _cache = {}

    def __init__(self, attributes={}):
        # Will be loaded as ndarray
        # TODO Jsonize here
        # Parse Ndarray serialised
        self._data = attributes.get("data", "{}")
        self._generation_id = int(attributes.get("orbit_generation_id", 0))





    def main(self):
        '''
        Main abstraction class for Input Manifold individual orbit processor
        :return:
        '''
        pass

    def generator(self):
        '''
        Main abstraction class for computation
        :return:
        '''
        pass

    def plot(self):
        '''
        Plotting specificity abstract class
        :return:
        '''
        pass

    @classmethod
    def clear_cache(cls):
        '''
        Method for cache control
        :return:
        '''
        cls._cache = {}
        cls._cache_remote_id = {}

    @property
    def data(self):
        if not self._data and self._generation_id:
            return None
        return self._data

    @data.setter
    def data(self, data):
        self._data = data

    def _insert_or_update_function_(self, is_ignore=False):
        '''
        Implement INSERT or UPDATE on DB
        :param is_ignore:
        :return:
        '''
        pass

    def save(self, is_ignore=False):
        '''
        Implement save
        :param is_ignore:
        :return:
        '''
        pass

    @classmethod
    def find(cls, id):

        # Retrieve processing unit from cache
        if id in cls._cache:
            return cls._cache[id]

        sql = f'''
            SELECT * FROM orbit_generation WHERE orbit_generation_id = {id}
        '''

        # TODO Connect to DB
        rows = []
        # rows = h.select(sql)

        if len(rows) == 0:
            return None

        return cls(rows[0])