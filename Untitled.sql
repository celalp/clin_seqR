CREATE TABLE "users" (
  "user_id" int,
  "username" int,
  "password" varchar,
  "created_at" datetime
);

CREATE TABLE "projects" (
  "project_id" int,
  "project_name" varchar,
  "description" varchar
);

CREATE TABLE "samples" (
  "sample_id" int,
  "sample_name" varchar,
  "description" int,
  "created_at" datetime
);

CREATE TABLE "sample_projects" (
  "id" int,
  "sample_id" int,
  "project_id" int
);

CREATE TABLE "users_projects" (
  "id" int,
  "project_id" int,
  "user_id" int
);

ALTER TABLE "sample_projects" ADD FOREIGN KEY ("sample_id") REFERENCES "samples" ("sample_id");

ALTER TABLE "sample_projects" ADD FOREIGN KEY ("project_id") REFERENCES "projects" ("project_id");

ALTER TABLE "users_projects" ADD FOREIGN KEY ("project_id") REFERENCES "projects" ("project_id");

ALTER TABLE "users_projects" ADD FOREIGN KEY ("user_id") REFERENCES "users" ("user_id");