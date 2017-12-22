# docker-compose-razi
Docker-compose to use compound database easily.

## Docker image
### db
PostgreSQL with [RDKit database cartridge](http://www.rdkit.org/docs/Cartridge.html) for compound database

- PostsgreSQL (version 9.5)
- RDKit database cartridge (latest)

### web
Jupyter Notebook for creating database table and using database with Python

- Python (version 3.6.0)
   - SQLAlchemy 
   - Pillow (version 4.3.0)
   - RDKit (version 2017.03.1)
   - [razi](https://github.com/rvianello/razi) (latest)
 
## How to use

### 1. git clone
```
git clone https://github.com/yamasakih/docker-compose-radi
cd docker-compose-django-radi
```

### 2. docker-compose up
```
docker-compose up -d
```

```
> $ docker-compose up -d
Creating network "dockercomposerazi_default" with the default driver
Creating dockercomposerazi_db_1 ... 
Creating dockercomposerazi_db_1 ... done
Creating dockercomposerazi_web_1 ... 
Creating dockercomposerazi_web_1 ... done
```

Please access `http://0.0.0.0:8888` or `http://localhost:8888`, then you can run Jupyter Notebook.

If you can't, please restart docker-compose.

```
> $ docker-compose restart
Restarting dockercomposedjangordkit_web_1 ... done
Restarting dockercomposedjangordkit_db_1  ... done
```

### 3. Creating the database
If you want, create database by command `createdb`.
You can use postgresql command by simply adding `docker-compose exec db` at the beginning.
(Change YOUR_DATABASE_NAME to created database name.)

```
> $ docker-compose exec db createdb -U postgres tutorial
```

### 4. Extend the database
Next, extend your database by command `create extension rdkit` for using RDKit database cartridge.
You can use postgresql command by simply adding `docker-compose exec db` at the beginning.
(Change YOUR_DATABASE_NAME to created database name or postgres (default database name).)

```
> $ docker-compose exec db psql -U postgres -c 'create extension rdkit' postgres
CREATE EXTENSION
```

### 5. Delete mount volume
This Docker-compose uses named volume. 
If it becomes unnecessary, delete it with the following command.

```
> $ docker volume list
DRIVER              VOLUME NAME
local               dockercomposerazi_dbdata
```

```
> $ docker volume rm dockercomposerazi_dbdata 
dockercomposerazi_dbdata
```

```
> $ docker volume list                        
DRIVER              VOLUME NAME
```

For more information, please see [tutorial](https://github.com/yamasakih/docker-compose-razi/tree/master/work/tutorial).
