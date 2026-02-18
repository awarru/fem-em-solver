#!/bin/bash
docker-compose -f docker/docker-compose.yml up -d
docker exec fem-em-solver python /workspace/examples/magnetostatics/01_straight_wire.py
