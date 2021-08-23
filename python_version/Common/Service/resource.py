import atexit
import hashlib
import os
import time
import urllib.request

import boto3
import magic
from boto3.s3.transfer import TransferConfig, create_transfer_manager

from Common.Base.singleton import SingletonMeta
from Common.Config.configurator import Configurator
from Common.Service.benchmark import Benchmark as b


class ResourceService(metaclass=SingletonMeta):
    temporal_images = []

    S3_KEY_ID = ''
    S3_ACCESS_KEY = ''
    S3_URL = ''

    def __init__(self):
        MB = 1024 ** 2
        config = TransferConfig(multipart_threshold=200 * MB)
        self.s3 = create_transfer_manager(
            boto3.client('s3', aws_access_key_id=Configurator().get('aws.s3.access_key_id', ''),
                         aws_secret_access_key=Configurator().get('aws.s3.secret_access_key',
                                                                  '')),
            config
        )
        self.base_url = Configurator().get('aws.s3.url', ResourceService)
        self.bucket_name = Configurator().get('aws.s3.name', 'aws-resources')

    def upload_resource(self, local_filename: str, resource_hash: str, filename: str, path: str, url: str = None,
                        bucket: str = None, is_public: bool = False):

        if bucket is None:
            bucket = self.bucket_name
        extra_args = {
            'ContentType': magic.from_file(local_filename, mime=True),
            'ContentDisposition': f'attachment; filename="{filename}"',
            'Metadata': {
                'filename': filename,
                'url': url or '-'
            }
        }
        if is_public:
            extra_args['ACL'] = 'public-read'
        return self.s3.upload(local_filename, bucket, f'{path}/{resource_hash}', extra_args=extra_args)

    def upload_resource_content(self, content, filename: str, path: str, bucket: str = None, is_public: bool = False):
        concurrent_filename = f'{filename}_{time.time_ns().real}'.encode()
        resource_hash = hashlib.md5(concurrent_filename).hexdigest()
        tmp_filename = os.getcwd() + f'/.cache/{resource_hash}'
        write_mode = 'wb+'
        if type(content) is str:
             write_mode = 'w+'
        f = open(tmp_filename, write_mode)
        f.write(content)
        f.close()
        self.temporal_images.append(tmp_filename)
        self.upload_resource(tmp_filename, resource_hash, filename, path, bucket=bucket, is_public=is_public)
        return resource_hash

    def upload_resource_content_safely(self, content, filename: str, path: str, bucket: str = None, is_public: bool = False, i=0, is_hashed=False):
        b.stat(f'Start uploading file at "{path}" at attempt #{i}', True)
        concurrent_filename = f'{filename}_{time.time_ns().real}'.encode()
        if not is_hashed:
            resource_hash = hashlib.md5(concurrent_filename).hexdigest()
        else:
            resource_hash = filename
        tmp_filename = os.getcwd() + f'/.cache/{resource_hash}'
        write_mode = 'wb+'
        if type(content) is str:
             write_mode = 'w+'
        f = open(tmp_filename, write_mode)
        f.write(content)
        f.close()
        self.temporal_images.append(tmp_filename)
        if bucket is None:
            bucket = self.bucket_name
        extra_args = {
            'ContentType': magic.from_file(tmp_filename, mime=True),
            'ContentDisposition': f'attachment; filename="{filename}"',
            'Metadata': {
                'filename': filename,
                'url': bucket or '-'
            }
        }
        if is_public:
            extra_args['ACL'] = 'public-read'
        aws_upload = self.s3.upload(tmp_filename, bucket, f'{path}/{resource_hash}', extra_args=extra_args)

        # wait for result
        result = aws_upload.result()

        # Check for download
        result_filename = self.download_resource(resource_hash=resource_hash, path=path, bucket=bucket)

        # Try again if nonexistant
        if result_filename is None:
            b.stat(f'Retry uploading file at path "{path}" at attempt #{i+1} as it previously failed', True)
            return self.upload_resource_content_safely(content=content, filename=filename, path=path, bucket=bucket, is_public=is_public, i=i+1)
        b.stat(f'Success uploading file at "{path}" with hash: {resource_hash} and name {filename}', True)
        return resource_hash

    def download_resource(self, resource_hash: str, path: str, bucket: str = None, real_name: str = None) -> str:
        filename = os.getcwd() + f'/.cache/{resource_hash if real_name is None else real_name}'
        skip_wait = False
        if bucket is None:
            bucket = self.bucket_name
        try:
            self.s3.download(bucket, f'{path}/{resource_hash}', filename)
        except BaseException as e:
            print(e)
            filename = None
            skip_wait = True
        finally:
            if filename is not None:
                self.temporal_images.append(filename)
        if not skip_wait:
            while not os.path.exists(filename):
                time.sleep(0.01)
        return filename

    def upload_image_from_url(self, resource_url: str, image_hash: str, bucket: str = None,
                              is_public: bool = False) -> bool:
        filename = os.getcwd() + f'/.cache/{image_hash}'
        name = resource_url.split('/')[len(resource_url.split('/')) - 1]
        try:
            urllib.request.urlretrieve(resource_url, filename)
            self.upload_resource(filename, image_hash, name, 'images', url=resource_url, bucket=bucket,
                                 is_public=is_public)
            success = True
        except BaseException as e:
            print(resource_url)
            print(e)
            success = False
        finally:
            self.temporal_images.append(filename)

        return success

    def get_image_from_hash(self, image_hash, bucket=None, path='images', real_name=None):
        return self.download_resource(image_hash, path, bucket=bucket, real_name=real_name)

    def clean_temporal_images(self, response=None):
        if len(self.temporal_images):
            print(f'\t[ResourceService] Cleaning up temporal images ({len(self.temporal_images)})')
            for filename in self.temporal_images:
                try:
                    if os.path.exists(filename):
                        os.remove(filename)
                except BaseException as e:
                    print(e)
        return response

    @staticmethod
    def clean_temporal_data():
        ResourceService().clean_temporal_images()


atexit.register(ResourceService.clean_temporal_data)


