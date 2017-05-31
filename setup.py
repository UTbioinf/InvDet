from skbuild import setup
import subprocess

setup(
    name="invdet",
    version=("0.2."+subprocess.check_output(["git", "rev-list", "--count", "f3e1a291253b5ac7ef2b1854adb879d115704d96..HEAD"])).rstrip(),
    description="Invertion Detection package",
    author="zijuexiansheng",
    license="MIT",
    packages=["invdet", "invdet.util"],
    package_dir={"invdet": "source/pyinvdet",
                 "invdet.util": "source/pyinvdet/util"},
    entry_points={
            'console_scripts': ['invdet=invdet.invdet_main:main']
        }
)
